import numpy as np
from astropy.convolution import convolve as apyconvolve
# from varinit import outvar
# from interplib_2 import get_zact

C_P = 1004
R_D = 287
R_V = 461 #Really 461.5, but this is what RAMS uses
EPS = R_D/R_V

def pvbuck(temp):
    tempc = temp-273.15; icetempc = np.where(tempc<0, tempc, np.nan)
    pveq_liq = 6.1121*np.exp((18.678-tempc/234.5)*(tempc/(257.14+tempc)))
    pveq_ice = 6.1115*np.exp((23.036-icetempc/333.7)*(icetempc/(279.82+icetempc)))
    return{"EqVaporPressLiq": pveq_liq, "EqVaporPressIce": pveq_ice}

def rams_pveq_ice(icetempc):
    c0= .6114327e+03; c1= .5027041e+02; c2= .1875982e+01
    c3= .4158303e-01; c4= .5992408e-03; c5= .5743775e-05
    c6= .3566847e-07; c7= .1306802e-09; c8= .2152144e-12

    x = np.where(icetempc>-80, icetempc, -80)
    x = np.where(-1*(np.isnan(icetempc))+1, x, np.nan)
    # x = np.where(tempc<0, tempc, 0)
    esif = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
    esif = np.where(x<0, esif, np.nan)
    return esif/100 #Convert from Pa to hPa
# In[10]:
def hydrovarsfunc(hydropath, refdens, hydromass, hydronum):
    import pandas as pd
    hydrodf = pd.read_csv(hydropath, delimiter = ",", index_col = 0);
    #Can I speed this up with multiprocessing? Each hydrometeor is independent, so this is embarassingly parallel
    hydroname = hydromass.name
    if hydroname == "Pris" or hydroname == "Snow":
        dfkey = f"{hydroname.lower()} col" #Get cfmas and pwmas properties from the hydrometeor table
    else:
        dfkey = hydroname.lower()
    numconc = hydronum*refdens #Convert number/kg to number/m^3
    mixvol = 1000*hydromass*refdens #Convert kg/kg to g/m^3
    vertint = (hydromass*refdens).integrate(coord = "z") #mm water equivalent of the given hydrometeor species
    diam = 10**6*(hydromass/(hydronum*hydrodf.loc[dfkey, "cfmas"]))**(1/hydrodf.loc[dfkey, "pwmas"]) #Hydrometeor mean diameter (microns)
    return {"Num": numconc, "MixVol": mixvol, "VertInt": vertint, "Diam": diam}

def get_rollvar(ywindow, xwindow, var3d, kerntype = "trikernel"):
    apxwindow = xwindow
    apywindow = ywindow
    if not ywindow%2:
        apywindow = ywindow+1 #For some reason, astropy needs odd-length windows
    if not xwindow%2:
        apxwindow = xwindow+1
    xwin1d = np.linspace(-apxwindow//2, apxwindow//2, apxwindow); ywin1d = np.linspace(-apxwindow//2, apywindow//2, apywindow)
    xwin2d, ywin2d = np.meshgrid(xwin1d, ywin1d, indexing = "xy")
    #Below are a few different kernel options for convolution.
    #Flat kernel weights everything within the window equally - I wouldn't recommend this, since it may cause artificial oscillations in the smoothed field
    flatkernel = np.ones((apywindow, apxwindow))
    #Trikernel creates a conical kernel, weighting cells near the center of the window more than cells toward the edge. This is a fairly well-behaved window, without much ringing
    trikernel = (1-(xwin2d**2+ywin2d**2)**(1/2)/(max((apxwindow-1)/2, (apywindow-1)/2)))
    trikernel = np.where(trikernel>=0, trikernel, 0)
    #Domekernel uses a parabolic convolutional kernel, with cells at the center given weight 1 and cells at the edge given weight 0. This has similar properties to trikernel, but a bit less smooth
    domekernel = (1-(xwin2d/((apxwindow-1)/2))**2-(ywin2d/((apywindow-1)/2))**2)
    domekernel = np.where(domekernel>=0, domekernel, 0)
    #Hornkernel uses a spire-like kernel that very heavily weights the cells near the center of the window, with almost no weight in the outher 2/3 of the window. I wouldn't recommend this, as it seems to cause oscillations similar to the boxcar window
    hornkernel = np.exp(-((xwin2d**2+ywin2d**2)**(1/2)/(max((apxwindow-1)/2, (apywindow-1)/2))))**4
    #Gausskernel is a gaussian kernel with weight of 1 in the center and lower weights on the edge. Unlike the trikernel and domekernel, this doesn't have a weight of 0 anywhere in the domain. gausskernel behaves similarly to trikernel and domekernel, but is slightly smoother
    gausskernel = np.exp(-(xwin2d**2+ywin2d**2)/(0.5*max((apxwindow-1)/2, (apywindow-1)/2)**2))
    kerndict = {"flatkernel": flatkernel, "trikernel": trikernel, "domekernel": domekernel, "hornkernel": hornkernel, "gausskernel": gausskernel}
    rollvar = np.zeros(var3d.shape)
    for zlev in range(rollvar.shape[0]):
        #Change the kernel below if you want to use something besides a trikernel
        rollvar[zlev,:,:] = apyconvolve(var3d[zlev,:,:], kernel = kerndict[kerntype], boundary = "extend", normalize_kernel = True, nan_treatment = "interpolate", preserve_nan = True)
    return rollvar


#In all of these routines, the "d" on the end of some variable names refers to the *data* for that variable. For example, "srftempd" is just the raw numpy array of terrain-following surface temperature. By contrast, "srftemp" is the outvar object which will go in the vardict array, containing various attributes like variable name, units, and information on the variables used to calculate it.
def get_srftemp(rawfile):
    #Return terrain-following surface temperature
    srftempd = rawfile["THETA"][1,1:-1,1:-1].values*rawfile["PI"][1,1:-1,1:-1].values/1004
    return srftempd

def get_srfpres(rawfile):
    #Return terrain-following surface pressure
    srfpresd = 1000*(rawfile["PI"][1,1:-1,1:-1].values/C_P)**(C_P/R_D)
    return srfpresd

def get_mslp(rawfile):
    #Return approximate sea-level pressure at the surface
    try:
        srftempd
    except:
        srftempd = rawfile["THETA"][1,1:-1,1:-1].values*rawfile["PI"][1,1:-1,1:-1].values/1004
    mslpd =  1000*(rawfile["PI"][1,1:-1,1:-1].values/1004)**(1004/287)*(srftempd/(srftempd+0.0065*rawfile["TOPT"][1:-1,1:-1].values))**(-9.80/(0.0065*287)) #Note the 0.0065 here is a 6.5 C/km lapse rate, which is also used by REVU's MSLP calculation. It may be worth revisiting this with a different lapse rate formulation in future.
    return mslpd

def get_mslt(rawfile):
    #Return approximate sea-level temperature at the surface (useful for figuring out which areas are warmer or colder adjusting for topography)
    try:
        srftempd
    except:
        srftempd = rawfile["THETA"][1,1:-1,1:-1].values*rawfile["PI"][1,1:-1,1:-1].values/1004
    msltd = srftempd+0.0065*rawfile["TOPT"][1:-1,1:-1].values #Again, using a 6.5 C/km lapse rate to approximate sea level temperature. Need a better approach than a constant lapse rate, though. Testing in my lake-effect sims shows that 9.8 C/km actually works much better, so it might depend on the situation
    return msltd

def get_condmix(vardict):
    #Return total condensate mixing ratio
    condmixd = vardict["RTP"].data-vardict["RV"].data
    return condmixd

def get_icemix(vardict):
    #Return total ice mixing ratio
    try:
        prismixd = vardict["RPP"].data
    except:
        prismixd = 0
    try:
        snowmixd = vardict["RSP"].data
    except:
        snowmixd = 0
    try:
        aggmixd = vardict["RAP"].data
    except:
        aggmixd = 0
    try:
        graupelmixd = vardict["RGP"].data
    except:
        graupelmixd = 0
    try:
        hailmixd = vardict["RHP"].data
    except:
        hailmixd = 0
    icemixd = prismixd+snowmixd+aggmixd+graupelmixd+hailmixd
    return icemixd

def get_liqmix(vardict):
    #Return total liquid mixing ratio
    try:
        cloudmixd = vardict["RCP"].data
    except:
        cloudmixd = 0
    try:
        drizzlemixd = vardict["rdP"].data
    except:
        drizzlemixd = 0
    try:
        rainmixd = vardict["RRP"].data
    except:
        rainmixd = 0
    liqmixd = cloudmixd+drizzlemixd+rainmixd
    return liqmixd

def get_thetav(vardict):
    #Return virtual potential temperature
    thetad = vardict["THETA"].data; vapormixd = vardict["RV"].data
    thetavd = thetad*(1+1/EPS*vapormixd)/(1+vapormixd)
    return thetavd

def get_thetarho(vardict):
    #Return density potential temperature
    try:
        thetad #save time/space if we already got this before
    except:
        thetad = vardict["THETA"].data
    try:
        vapormixd
    except:
        vapormixd = vardict["RV"].data; 
    watermixd = vardict["RTP"].data
    thetarhod = thetad*(1+vapormixd/EPS)/(1+watermixd)
    return thetarhod

def get_pressure(vardict):
    #Return pressure field
    pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
    return pressured

def get_temp(vardict):
    #Return temperature field
    try:
        thetad
    except:
        thetad = vardict["THETA"].data
    tempd = thetad*vardict["PI"].data/C_P
    return tempd

def get_tempv(vardict):
    #Return virtual temperature
    try:
        thetavd
    except:
        thetad = vardict["THETA"].data; vapormixd = vardict["RV"].data
        thetavd = thetad*(1+1/EPS*vapormixd)/(1+vapormixd)
    tempvd = thetavd*vardict["PI"].data/C_P
    return tempvd

def get_temprho(vardict):
    #Return density temperature
    try:
        thetarhod
    except:
        thetarhod = vardict["THETA"].data*(1+vardict["RV"].data/EPS)/(1+vardict["RTP"].data)
    temprhod = thetarhod*vardict["PI"].data/C_P
    return temprhod

def get_rho(vardict):
    #Return air density (using ideal gas law, doesn't account for condensate loading)
    try:
        pressured
    except:
        pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
    try:
        tempvd
    except:
        try:
            thetavd
        except:
            try:
                thetad
            except:
                thetad = vardict["THETA"].data; vapormixd = vardict["RV"].data
            thetavd = thetad*(1+1/EPS*vapormixd)/(1+vapormixd)
        tempvd = thetavd*vardict["PI"].data/C_P
    rhod = 100*pressured/(tempvd*R_D) #Factor of 100 here is because pressured is stored in hPa, while rhod is in SI units (kg/m^3), so need to convert hPa to Pa
    return rhod

def get_vaporpres(vardict):
    #Return vapor pressure
    try:
        pressured
    except:
        pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
    vaporpresd = vardict["RV"].data*pressured/(EPS+vardict["RV"].data)
    return vaporpresd

def get_pwat(rawfile, rawzsub):
    #Return precipitable water in each column
    pressureraw = 1000*(rawfile["PI"][1:-1,1:-1,1:-1].values/C_P)**(C_P/R_D)
    vaporpresraw = rawfile["RV"][1:-1,1:-1,1:-1].values*pressureraw/(EPS+rawfile["RV"][1:-1,1:-1,1:-1].values)
    tempraw = rawfile["THETA"][1:-1,1:-1,1:-1].values*rawfile["PI"][1:-1,1:-1,1:-1].values/C_P
    pwatd =  np.trapz(rawfile["RV"][1:-1,1:-1,1:-1].values*100*(pressureraw-vaporpresraw)/(tempraw*R_D), rawzsub, axis = 0)
    return pwatd
    #Note that even if we're using cartesian or pressure interpolation, or a shorter subset of the sigma-z domain, we're using the full sigma-z grid for this calculation. This is because we want to capture all the water vapor up to the model top, which might be above the domain subset the user has chosen.

def get_dewpoint(vardict):
    try:
        vaporpresd
    except:
        try:
            pressured
        except:
            pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
        vaporpresd = vardict["RV"].data*pressured/(EPS+vardict["RV"].data)
    dewpointd = (4717.2-35.86*np.log(vaporpresd/6.1078))/(17.269-np.log(vaporpresd/6.1078))
    return dewpointd

def get_cloudtopheight(rawfile, rawzsub):
    cloudtopheightd = np.nanmax(np.where(rawfile["RCP"][1:-1,1:-1,1:-1].values+rawfile["RPP"][1:-1,1:-1,1:-1].values>1e-4, rawzsub, np.nan), axis = 0)
    #Again, like above, using the full sigma-z grid for the same reason as above - the user's analyais domain top might be below the heights of the highest clouds
    return cloudtopheightd

def get_cloudbaseheight(rawfile, rawzsub):
    cloudbaseheightd = np.nanmin(np.where(rawfile["RCP"][1:-1,1:-1,1:-1].values+rawfile["RPP"][1:-1,1:-1,1:-1].values>1e-4, rawzsub, np.nan), axis = 0)
    return cloudbaseheightd

def get_cloudtoppres(rawfile):
    try:
        pressureraw
    except:
        pressureraw = 1000*(rawfile["PI"][1:-1,1:-1,1:-1].values/C_P)**(C_P/R_D)
    cloudtoppresd = np.nanmin(np.where(rawfile["RCP"][1:-1,1:-1,1:-1].values+rawfile["RPP"][1:-1,1:-1,1:-1].values>1e-4, pressureraw, np.nan), axis = 0)
    return cloudtoppresd

def get_cloudbasepres(rawfile):
    try:
        pressureraw
    except:
        pressureraw = 1000*(rawfile["PI"][1:-1,1:-1,1:-1].values/C_P)**(C_P/R_D)
    cloudbasepresd = np.nanmax(np.where(rawfile["RCP"][1:-1,1:-1,1:-1].values+rawfile["RPP"][1:-1,1:-1,1:-1].values>1e-4, pressureraw, np.nan), axis = 0)
    return cloudbasepresd

def get_wadv_h(vardict, gridprops):
    wadvhd = -vardict["VC"].data*np.gradient(vardict["WC"].data, gridprops["dx"], axis = 1)-vardict["UC"].data*np.gradient(vardict["WC"].data, gridprops["dx"], axis = 2)
    return wadvhd

def get_wadv_v(vardict, ccoords):
    wadvvd = -vardict["WC"].data*np.gradient(vardict["WC"].data, ccoords["z"], axis = 0)
    return wadvvd

def get_horizdiv(vardict, gridprops):
    divd = np.gradient(vardict["UC"].data, gridprops["dx"], axis = 2)+np.gradient(vardict["VC"].data, gridprops["dx"], axis = 1)
    return divd

def get_xvort(vardict, gridprops, ccoords):
    xvortd = np.gradient(vardict["WC"].data, gridprops["dx"], axis = 1)-np.gradient(vardict["VC"].data, ccoords["z"].values, axis = 0)
    return xvortd

def get_yvort(vardict, gridprops, ccoords):
    yvortd = np.gradient(vardict["UC"].data, ccoords["z"].values, axis = 0)-np.gradient(vardict["WC"].data, gridprops["dx"], axis = 2)
    return yvortd

def get_zvortrel(vardict, gridprops):
    zvortreld = np.gradient(vardict["VC"].data, gridprops["dx"], axis = 2) - np.gradient(vardict["UC"].data, gridprops["dx"], axis = 1)
    return zvortreld

def get_zvortabs(vardict, gridprops, ccoords):
    zvortabsd = np.gradient(vardict["VC"].data, gridprops["dx"], axis = 2) - np.gradient(vardict["UC"].data, gridprops["dx"], axis = 1) + 2*7.292*10**(-5)*np.sin(ccoords["lat2d"].values[None,:,:]*np.pi/180)
    return zvortabsd
#The [None,:,:] on lat2d extends it to a phony axis in the z dimension so that numpy can broadcast it to 3D and match the shape of UC and VC

def get_liqsatfrac(vardict):
    #Returns the fraction of equilibrium vapor pressure with respect to water of the current vapor pressure (0 = 0% humidity, 1 = 100% humidity, 2 = 200% humidity, etc.)
    try:
        vaporpresd
    except:
        try:
            pressured
        except:
            pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
        vaporpresd = vardict["RV"].data*pressured/(EPS+vardict["RV"].data)
    try:
        tempd
    except:
        tempd = vardict["THETA"].data*vardict["PI"].data/C_P
    tempc = tempd-273.15
    pveq_liq = 6.1121*np.exp((18.678-tempc/234.5)*(tempc/(257.14+tempc)))
    liqsatfracd = vaporpresd/pveq_liq
    return liqsatfracd

def get_liqvapordeficit(vardict):
    #Returns the difference between the current vapor pressure and the equilibrium vapor pressure with respect to water at each gridpoint.
    try:
        liqsatfracd
    except:
        try:
            vaporpresd
        except:
            try:
                pressured
            except:
                pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
            vaporpresd = vardict["RV"].data*pressured/(EPS+vardict["RV"].data)
        try:
            tempd
        except:
            tempd = vardict["THETA"].data*vardict["PI"].data/C_P
    tempc = tempd-273.15
    pveq_liq = 6.1121*np.exp((18.678-tempc/234.5)*(tempc/(257.14+tempc)))
    liqsatfracd = vaporpresd/pveq_liq
    liqvapordeficitd = pveq_liq*(1-liqsatfracd)
    return liqvapordeficitd

def get_icesatfrac(vardict):
    try:
        vaporpresd
    except:
        try:
            pressured
        except:
            pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
        vaporpresd = vardict["RV"].data*pressured/(EPS+vardict["RV"].data)
    try:
        tempd
    except:
        tempd = vardict["THETA"].data*vardict["PI"].data/C_P
    tempc = tempd-273.15
    icetempc = np.where(tempc<0, tempc, np.nan)
    pveq_ice = rams_pveq_ice(icetempc)
    icesatfracd = vaporpresd/pveq_ice
    return icesatfracd

def get_icevapordeficit(vardict):
    try:
        icesatfracd
    except:
        try:
            vaporpresd
        except:
            try:
                pressured
            except:
                pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
            vaporpresd = vardict["RV"].data*pressured/(EPS+vardict["RV"].data)
        try:
            tempd
        except:
            tempd = vardict["THETA"].data*vardict["PI"].data/C_P
        tempc = tempd-273.15
        icetempc = np.where(tempc<0, tempc, np.nan)
        pveq_ice = rams_pveq_ice(icetempc)
        icesatfracd = vaporpresd/pveq_ice
    icevapordeficitd = pveq_ice*(1-icesatfracd)
    return icevapordeficitd

def get_superliqmix(vardict):
    #Return the difference between the current vapor mixing ratio and the saturation vapor mixing ratio with respect to liquid at each gridpoint
    try:
        pveq_liq
    except:
        try:
            vaporpresd
        except:
            try:
                pressured
            except:
                pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
            vaporpresd = vardict["RV"].data*pressured/(EPS+vardict["RV"].data)
        try:
            tempd
        except:
            tempd = vardict["THETA"].data*vardict["PI"].data/C_P
        tempc = tempd-273.15
        pveq_liq = 6.1121*np.exp((18.678-tempc/234.5)*(tempc/(257.14+tempc)))
    try:
        rhodd #rhodd is the *dry air density data*, not the total air density data (rhod). This is because mixing ratio is vapor density/dry air density, not vapor density/total density
    except:
        rhodd = (pressured-vaporpresd)/(R_D*tempd)
    superliqmixd = 100*(vaporpresd-pveq_liq)/(R_V*tempd*rhodd) #This divides the density of the excess vapor (in kg/m^3, and then divides by the dry air density to get a mixing ratio)
    return superliqmixd

def get_supericemix(vardict):
    try:
        pveq_ice
    except:
        try:
            vaporpresd
        except:
            try:
                pressured
            except:
                pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
            vaporpresd = vardict["RV"].data*pressured/(EPS+vardict["RV"].data)
        try:
            tempd
        except:
            tempd = vardict["THETA"].data*vardict["PI"].data/C_P
        tempc = tempd-273.15
        icetempc = np.where(tempc<0, tempc, np.nan)
        pveq_ice = rams_pveq_ice(icetempc)
    try:
        rhodd
    except:
        rhodd = (pressured-vaporpresd)/(R_D*tempd)
    supericemixd = 100*(vaporpresd-pveq_ice)/(R_V*tempd*rhodd)
    return supericemixd

def get_cloudnum(vardict):
    cloudnumd = vardict["CCP"].data*vardict["DN0"].data
    return cloudnumd

def get_drizzlenum(vardict):
    drizzlenumd = vardict["CDP"].data*vardict["DN0"].data
    return drizzlenumd

def get_rainnum(vardict):
    rainnumd = vardict["CRP"].data*vardict["DN0"].data
    return rainnumd

def get_prisnum(vardict):
    prisnumd =  vardict["CPP"].data*vardict["DN0"].data
    return prisnumd

def get_snownum(vardict):
    snownumd = vardict["CSP"].data*vardict["DN0"].data
    return snownumd

def get_aggregatenum(vardict):
    aggnumd = vardict["CAP"].data*vardict["DN0"].data
    return aggnumd

def get_graupelnum(vardict):
    graupelnumd = vardict["CGP"].data*vardict["DN0"].data
    return graupelnumd

def get_hailnum(vardict):
    hailnumd = vardict["CHP"].data*vardict["DN0"].data
    return hailnumd

def get_platenum(vardict):
    platenumd = vardict["CIPP"].data*vardict["DN0"].data
    return platenumd

def get_columnnum(vardict):
    columnnumd = vardict["CICP"].data*vardict["DN0"].data
    return columnnumd

def get_dendritenum(vardict):
    dendritenumd = vardict["CIDP"].data*vardict["DN0"].data
    return dendritenumd

def get_cloudwatercontent(vardict):
    cloudcond = vardict["RCP"].data*vardict["DN0"].data
    return cloudcond

def get_drizzlewatercontent(vardict):
    drizcond = vardict["RDP"].data*vardict["DN0"].data
    return drizcond

def get_rainwatercontent(vardict):
    raincond = vardict["RRP"].data*vardict["DN0"].data
    return raincond

def get_priswatercontent(vardict):
    priscond = vardict["RPP"].data*vardict["DN0"].data
    return priscond

def get_snowwatercontent(vardict):
    snowcond = vardict["RSP"].data*vardict["DN0"].data
    return snowcond

def get_aggregatewatercontent(vardict):
    aggcond = vardict["RAP"].data*vardict["DN0"].data
    return aggcond

def get_graupelwatercontent(vardict):
    graupelcond = vardict["RGP"].data*vardict["DN0"].data
    return graupelcond

def get_hailwatercontent(vardict):
    hailcond = vardict["RHP"].data*vardict["DN0"].data
    return hailcond

def get_platewatercontent(vardict):
    platecond = vardict["RIPP"].data*vardict["DN0"].data
    return platecond

def get_columnwatercontent(vardict):
    columncond = vardict["RICP"].data*vardict["DN0"].data
    return columncond

def get_dendritewatercontent(vardict):
    dendritecond = vardict["RIDP"].data*vardict["DN0"].data
    return dendritecond

def get_cloudvertint(rawfile, rawzsub):
    cloudvertintd = np.trapz(rawfile["RCP"][1:-1,1:-1,1:-1].values*rawfile["DN0"][1:-1,1:-1,1:-1].values, rawzsub, axis = 0)
    return cloudvertintd

def get_drizzlevertint(rawfile, rawzsub):
    drizvertintd = np.trapz(rawfile["RDP"][1:-1,1:-1,1:-1].values*rawfile["DN0"][1:-1,1:-1,1:-1].values, rawzsub, axis = 0)
    return drizvertintd

def get_rainvertint(rawfile, rawzsub):
    rainvertintd = np.trapz(rawfile["RRP"][1:-1,1:-1,1:-1].values*rawfile["DN0"][1:-1,1:-1,1:-1].values, rawzsub, axis = 0)
    return rainvertintd

def get_prisvertint(rawfile, rawzsub):
    prisvertintd = np.trapz(rawfile["RPP"][1:-1,1:-1,1:-1].values*rawfile["DN0"][1:-1,1:-1,1:-1].values, rawzsub, axis = 0)
    return prisvertintd

def get_snowvertint(rawfile, rawzsub):
    snowvertintd = np.trapz(rawfile["RSP"][1:-1,1:-1,1:-1].values*rawfile["DN0"][1:-1,1:-1,1:-1].values, rawzsub, axis = 0)
    return snowvertintd

def get_aggregatevertint(rawfile, rawzsub):
    aggvertintd = np.trapz(rawfile["RAP"][1:-1,1:-1,1:-1].values*rawfile["DN0"][1:-1,1:-1,1:-1].values, rawzsub, axis = 0)
    return aggvertintd

def get_graupelvertint(rawfile, rawzsub):
    graupelvertintd = np.trapz(rawfile["RGP"][1:-1,1:-1,1:-1].values*rawfile["DN0"][1:-1,1:-1,1:-1].values, rawzsub, axis = 0)
    return graupelvertintd

def get_hailvertint(rawfile, rawzsub):
    hailvertintd = np.trapz(rawfile["RHP"][1:-1,1:-1,1:-1].values*rawfile["DN0"][1:-1,1:-1,1:-1].values, rawzsub, axis = 0)
    return hailvertintd

def get_platevertint(rawfile, rawzsub):
    platevertintd = np.trapz(rawfile["RIPP"][1:-1,1:-1,1:-1].values*rawfile["DN0"][1:-1,1:-1,1:-1].values, rawzsub, axis = 0)
    return platevertintd

def get_columnvertint(rawfile, rawzsub):
    columnvertintd = np.trapz(rawfile["RICP"][1:-1,1:-1,1:-1].values*rawfile["DN0"][1:-1,1:-1,1:-1].values, rawzsub, axis = 0)
    return columnvertintd

def get_dendritevertint(rawfile, rawzsub):
    dendritevertintd = np.trapz(rawfile["RIDP"][1:-1,1:-1,1:-1].values*rawfile["DN0"][1:-1,1:-1,1:-1].values, rawzsub, axis = 0)
    return dendritevertintd

def get_liqvertint(rawfile, rawzsub):
    try:
        liqmixd
    except:
        try:
            cloudmixd = rawfile["RCP"][1:-1,1:-1,1:-1].values
        except:
            cloudmixd = 0
        try:
            drizzlemixd = rawfile["RDP"][1:-1,1:-1,1:-1].values
        except:
            drizzlemixd = 0
        try:
            rainmixd = rawfile["RRP"][1:-1, 1:-1, 1:-1].values
        except:
            rainmixd = 0
        liqmixd = (cloudmixd+drizzlemixd+rainmixd)
    liqvertintd = np.trapz(liqmixd*(rawfile["DN0"][1:-1,1:-1,1:-1].values), rawzsub, axis = 0)
    return liqvertintd

def get_icevertint(rawfile, rawzsub):
    try:
        icemixd
    except:
        try:
            prismixd = rawfile["RPP"][1:-1,1:-1,1:-1].values
        except:
            prismixd = 0
        try:
            snowmixd = rawfile["RSP"][1:-1,1:-1,1:-1].values
        except:
            snowmixd = 0
        try:
            aggmixd = rawfile["RAP"][1:-1,1:-1,1:-1].values
        except:
            aggmixd = 0
        try:
            graupelmixd = rawfile["RGP"][1:-1,1:-1,1:-1].values
        except:
            graupelmixd = 0
        try:
            hailmixd = rawfile["RHP"][1:-1,1:-1,1:-1].values
        except:
            hailmixd = 0
        icemixd = (prismixd+snowmixd+aggmixd+graupelmixd+hailmixd)
    icevertintd = np.trapz(icemixd*(rawfile["DN0"][1:-1,1:-1,1:-1].values), rawzsub, axis = 0)
    return icevertintd

def get_condvertint(rawfile, rawzsub):
    try:
        condmixd
    except:
        condmixd = rawfile["RTP"][1:-1,1:-1,1:-1].values - rawfile["RV"][1:-1,1:-1,1:-1].values
    condvertintd = np.trapz(condmixd*(rawfile["DN0"][1:-1,1:-1,1:-1].values), rawzsub, axis = 0)
    return condvertintd

def get_clouddiam(vardict, hydrodf):
    clouddiamd = 1000*(vardict["RCP"].data/(vardict["CCP"].data*hydrodf.loc["cloud", "cfmas"]))**(1/hydrodf.loc["cloud", "pwmas"])
    return clouddiamd

def get_drizzlediam(vardict, hydrodf):
    drizdiamd = 1000*(vardict["RDP"].data/(vardict["CDP"].data*hydrodf.loc["drizzle", "cfmas"]))**(1/hydrodf.loc["drizzle", "pwmas"])
    return drizdiamd

def get_raindiam(vardict, hydrodf):
    raindiamd = 1000*(vardict["RRP"].data/(vardict["CRP"].data*hydrodf.loc["rain", "cfmas"]))**(1/hydrodf.loc["rain", "pwmas"])
    return raindiamd

def get_prisdiam(vardict, hydrodf):
    prisdiamd = 1000*(vardict["RPP"].data/(vardict["CPP"].data*hydrodf.loc["pris col", "cfmas"]))**(1/hydrodf.loc["pris col", "pwmas"])
    return prisdiamd

def get_snowdiam(vardict, hydrodf):
    snowdiamd = 1000*(vardict["RSP"].data/(vardict["CSP"].data*hydrodf.loc["snow col", "cfmas"]))**(1/hydrodf.loc["snow col", "pwmas"])
    return snowdiamd

def get_aggregatediam(vardict, hydrodf):
    aggdiamd = 1000*(vardict["RAP"].data/(vardict["CAP"].data*hydrodf.loc["agg", "cfmas"]))**(1/hydrodf.loc["agg", "pwmas"])
    return aggdiamd

def get_graupeldiam(vardict, hydrodf):
    graupeldiamd = 1000*(vardict["RGP"].data/(vardict["CGP"].data*hydrodf.loc["graupel", "cfmas"]))**(1/hydrodf.loc["graupel", "pwmas"])
    return graupeldiamd

def get_haildiam(vardict, hydrodf):
    haildiamd = 1000*(vardict["RHP"].data/(vardict["CHP"].data*hydrodf.loc["hail", "cfmas"]))**(1/hydrodf.loc["hail", "pwmas"])
    return haildiamd

def get_rhoprime(vardict, window, kernname):
    ywindow = window["ywindow"]; xwindow = window["xwindow"]
    try:
        rhod
    except:
        try:
            pressured
        except:
            pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
        try:
            tempvd
        except:
            try:
                thetavd
            except:
                try:
                    thetad
                except:
                    thetad = vardict["THETA"].data; vapormixd = vardict["RV"].data
                thetavd = thetad*(1+1/EPS*vapormixd)/(1+vapormixd)
            tempvd = thetavd*vardict["PI"].data/C_P
        rhod = 100*pressured/(tempvd*R_D)
    try:
        condmixd
    except:
        condmixd = vardict["RTP"].data-vardict["RV"].data
    rollrho = get_rollvar(ywindow, xwindow, rhod, kerntype = kernname)
    rhoprimed = rhod*(1+condmixd)-rollrho
    return rhoprimed

def get_rhobuoy(vardict, window, kernname):
    ywindow = window["ywindow"]; xwindow = window["xwindow"]
    try:
        condmixd
    except:
        condmixd = vardict["RTP"].data-vardict["RV"].data
    try:
        rhoprimed
    except:
        try:
            rhod
        except:
            try:
                pressured
            except:
                pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
            try:
                tempvd
            except:
                try:
                    thetavd
                except:
                    try:
                        thetad
                    except:
                        thetad = vardict["THETA"].data; vapormixd = vardict["RV"].data
                    thetavd = thetad*(1+1/EPS*vapormixd)/(1+vapormixd)
                tempvd = thetavd*vardict["PI"].data/C_P
            rhod = 100*pressured/(tempvd*R_D)
        rollrho = get_rollvar(ywindow, xwindow, rhod, kerntype = kernname)
        rhoprimed = rhod*(1+condmixd)-rollrho
    rhobuoyd = -9.80*(rhoprimed/rhod)
    return rhobuoyd

def get_thetabuoy(vardict, window, kernname):
    ywindow = window["ywindow"]; xwindow = window["xwindow"]
    try:
        thetarhod
    except:
        thetarhod = vardict["THETA"].data*(1+1/EPS*vardict["RV"].data)/(1+vardict["RTP"].data)
    rollthetarho = get_rollvar(ywindow, xwindow, thetarhod, kerntype = kernname)
    thetarhobuoyd = 9.80*((thetarhod-rollthetarho)/rollthetarho)
    return thetarhobuoyd

def get_thermalbuoy(vardict, window, kernname):
    ywindow = window["ywindow"]; xwindow = window["xwindow"]
    try:
        thetad
    except:
        thetad = vardict["THETA"].data
    rolltheta = get_rollvar(ywindow, xwindow, thetad, kerntype = kernname)
    thermbuoyd = 9.8*(thetad-rolltheta)/thetad
    return thermbuoyd

def get_vaporbuoy(vardict, window, kernname):
    ywindow = window["ywindow"]; xwindow = window["xwindow"]
    rollvap = get_rollvar(ywindow, xwindow, vardict["RV"].data, kerntype = kernname)
    vapbuoyd = 9.8*0.6077*(vardict["RV"].data-rollvap)
    return vapbuoyd

def get_condbuoy(vardict):
    try:
        condmixd
    except:
        condmixd = vardict["RTP"].data-vardict["RV"].data
    condbuoyd = -9.8*condmixd
    return condbuoyd

def get_pprimebuoy(vardict, window, kernname):
    ywindow = window["ywindow"]; xwindow = window["xwindow"]
    try:
        pressured
    except:
        pressured = 1000*(vardict["PI"].data/C_P)**(C_P/R_D)
    rollpressure = get_rollvar(ywindow, xwindow, pressured, kerntype = kernname)
    ppbuoyd = -9.80*0.714*(pressured-rollpressure)/pressured
    return ppbuoyd

def get_vppgf(vardict, window, kernname, ccoords):
    ywindow = window["ywindow"]; xwindow = window["xwindow"]
    try:
        thetarhod
    except:
        thetarhod = vardict["THETA"].data*(1+vardict["RV"].data/EPS)/(1+vardict["RTP"].data)
    rollpi = get_rollvar(ywindow, xwindow, vardict["PI"].data, kerntype = kernname)
    vppgfd = -thetarhod*(np.gradient(vardict["PI"].data-rollpi, ccoords["z"].values, axis = 0))
    return vppgfd











