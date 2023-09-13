#!/usr/bin/env python
# coding: utf-8

# In[2]:


import xarray as xr
import numpy as np
import pandas as pd
from time import perf_counter
from astropy.convolution import convolve as apyconvolve #Turns out astropy has a really convenient function for doing rolling means with masked data
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from varinit import outvar
from interplib_2 import get_zact
from derivcalcs import *
from dvardict import derivvarinit

'''
List of available derived variables. 
Format VerboseName: RAMSNAME --> Description of the variablee
SrfTemp: TS --> Terrain-following temperature at the surface, in K
SrfPres: PS --> Terrain-following pressure at the surface, in hPA
MSLP: MSLP --> Approximate sea-level pressure at each x-y gridpoint. Found by a hydrostatic approximation using a 6.5K/km lapse rate
MSLT: MSLT --> Approximate sea-level temperature at each x-y gridpoint. Found by interpolating surface pressure down from terrain height to sea level using a 6.5K/km lapse rate
CondMix: RTC --> Total condensate mixing ratio, in g/kg
IceMix: RTI --> Total ice condensate mixing ratio, in g/kg
LiqMix: RTL --> Total liquid condensate mixing ratio, in g/kg
ThetaV: THV --> Virtual potential temperature, in K
ThetaRho: THR --> Density potential temperature, in K
Pressure: P --> Pressure, in hPa
Temperature: T --> Temperature, in K
TempV: TV --> Virtual temperature, in K
TempRho: TR --> Density temperature, in K
RHO: RHO --> Density, in kg/m^3
VaporPressure: E --> Water vapor pressure, in hPa
PWAT: PWAT --> Precipitable water in each x-y column, in mm
Dewpoint: TD --> Dewpoint temperature, in K
CloudTopHeight: ZCT --> Highest point in each x-y column with cloud water or pristine ice, in m
CloudBaseHeight: ZCB --> Lowest point in each x-y column with cloud water or pristine ice, in m
CloudTopPressure: PCT --> Pressure at highest point in each x-y column with cloud water or pristine ice, in hPa
CloudBasePressure: PCB --> Pressure at lowest point in each x-y column with cloud water or pristine ice, in hPa
WAdv_H: WADVH --> Horizontal advection of vertical velocity, in m/s^2
WAdv_V: WADVV --> Vertical advection of vertical velocity, in m/s^2
XVort: VOX --> Magnitude and sign of vorticity in the +X (eastward) direction, in s^-1
YVort: VOY --> Magnitude and sign of vorticity in the +Y (northward) direction, in s^-1
ZVortRel: VOZR --> Magnitude and sign of relative vorticity in the +Z (upward) direction, in s^-1
ZVortAbs: VOZA --> Magnitude and sign of absolute vorticity in the +Z (upward) direction, in s^-1
LiqSatFrac: SL --> Water Vapor Pressure Fraction of Equilibrium Vapor Pressure with Respect to Liquid (unitless)
IceSatFrac: SI --> Water Vapor Pressure Fraction of Equilibrium Vapor Pressure with Respect to Ice (unitless)
LiqVaporDeficit: DVL --> Difference Between Equilibrium Vapor Pressure With Respect to Liquid and Current Water Vapor Pressure, in hPa
IceVaporDeficit: DVI --> Difference Between Equilibrium Vapor Pressure with Respect to Ice and Current Water Vapor Pressure, in hPa
SuperLiqMix: RVSI --> Mixing Ratio Difference from Equilibrium Water Vapor Mixing Ratio with Respect to Liquid, in g/kg
SuperIceMix: RVSI --> Mixing Ratio Difference from Equilibrium Water Vapor Mixing Ratio with Respect to Ice, in g/kg
CloudNum: NCP --> Number concentration of cloud droplets, in #/m^3
DrizNum: NDP --> Number concentration of drizzle, in #/m^3
RainNum: NRP --> Number concentration of rain, in #/m^3
PrisNum: NPP --> Number concentration of pristine ice, in #/m^3
SnowNum: NSP --> Number concentration of snow, in #/m^3
AggNum: NAP --> Number concentration of aggregates, in #/m^3
GraupelNum: NGP --> Number concentration of graupel, in #/m^3
HailNum: NHP --> Number concentration of hail, in #/m^3
PlateNum: NIPP --> Number concentration of plate ice, in #/m^3
ColumnNum: NICP --> Number concentration of column ice, in #/m^3
DendriteNum: NIDP --> Number concentration of dendrite ice, in #/m^3
CloudWaterContent: CWC
DrizzleWaterContent: DWC
RainWaterContent: RWC
PrisWaterContent: PWC
SnowWaterContent: SWC
AggWaterContent: AWC
GraupelWaterContent: GWC
HailWaterContent: HWC,
PlateWaterContent: IPWC
ColumnWaterContent: ICWC
DendriteWaterContent: IDWC
VertIntCloud: VIC
VertIntDrizzle: VID
VertIntRain: VIR
VertIntPris: VIP
VertIntSnow: VIS
VertIntAgg: VIA
VertIntGraupel: VIG
VertIntHail: VIH
VertIntPlate: VIIP
VertIntColumn: VIIC
VertIntDendrite: VIID
VertIntLiq: VIL
VertIntIce: VII
CloudDiam: DCP
DrizzleDiam: DDP
RainDiam: DRP
PrisDiam: DPP
SnowDiam: DSP
AggDiam: DAP
GraupelDiam: DGP
HailDiam: DHP
RhoPrime: RHOP
RhoBuoy: BUOY_RHO
ThetaRhoBuoy: BUOY_THETA
ThermalBuoy: BUOY_TEMP
VapBuoy: BUOY_VAP
CondBuoy: BUOY_COND
PPrimeBuoy: BUOY_PPRIME
VPPGF: VPPGF
'''

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
def hydrovarsfunc(hydrodf, refdens, hydromass, hydronum):
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

def get_derivedvars(vardict, derivvars, rawfile, gridprops, window, kernname, rnameflag, hydropath, ccoords = None):
    rawzsub = get_zact(gridprops, rawfile["TOPT"])[1:-1,1:-1,1:-1] #This gives sigma-z coords for the entire domain of the raw file (still cutting off model top, bottom, and horizontal boundaries as intended). Used for vertically intregated quantities and cloud tops/bottoms
    varkeys = list(vardict.keys())
    # ywindow = window["ywindow"]; xwindow = window["xwindow"]
    hydrodf = pd.read_csv(hydropath, delimiter = ",", index_col = 0); #print(hydroparams)
    derivvardict = derivvarinit()
    print(f"Finding Derived Variables {derivvars}!")
    #This for loop adds all derived variable objects (not the data) to vardict
    for dvarname in derivvars:
        if all(invar in varkeys for invar in derivvardict[dvarname].invar) and dvarname not in ["RTL", "RTI", "VIL", "VII"]: #Special cases for liquid and ice mixing ratios and vertically integrated quantities, because not all hydrometeor categories may be present
            vardict[dvarname] = derivvardict[dvarname]
        elif any(invar in varkeys for invar in derivvardict[dvarname].invar) and dvarname in ["RTL", "RTI", "VIL", "VII"]:
            vardict[dvarname] = derivvardict[dvarname]
        elif dvarname not in ["RTL", "RTI", "VIL", "VII"] and rnameflag == 0:
            print(f"Not all necessary variables found for computation of {derivvardict[dvarname].varname}! {derivvardict[dvarname].varname} needs variables {derivvardict[dvarname].invar}!")
        elif dvarname not in ["RTL", "RTI", "VIL", "VII"] and rnameflag == 1:
            print(f"Not all necessary variables found for computation of {derivvardict[dvarname].ramsname}! {derivvardict[dvarname].ramsname} needs variables {derivvardict[dvarname].invar}!")
        elif dvarname in ["RTL", "RTI", "VIL", "VII"] and rnameflag == 0:
            print(f"Not all necessary variables found for computation of {derivvardict[dvarname].varname}! {derivvardict[dvarname].varname} needs at least one of the following variables: {derivvardict[dvarname].invar}!")
        elif dvarname in ["RTL", "RTI", "VIL", "VII"] and rnameflag == 1:
            print(f"Not all necessary variables found for computation of {derivvardict[dvarname].ramsname}! {derivvardict[dvarname].ramsname} needs at least one of the following variables: {derivvardict[dvarname].invar}!")
        
        

    #Now, actually calculate the data associated with each derived variable
    if "TS" in vardict.keys():
        srftempd = get_srftemp(rawfile)
        vardict["TS"].data = srftempd
        print(f"Found variable {[vardict['TS'].varname, vardict['TS'].ramsname][rnameflag]}!")

    if "PS" in vardict.keys():
        srfpresd = get_srfpres(rawfile)
        vardict["PS"].data = srfpresd
        print(f"Found variable {[vardict['PS'].varname, vardict['PS'].ramsname][rnameflag]}!")

    if "MSLP" in vardict.keys():
        mslpd = get_mslp(rawfile)
        vardict["MSLP"].data = mslpd
        print(f"Found variable {[vardict['MSLP'].varname, vardict['MSLP'].ramsname][rnameflag]}!")

    if "MSLT" in vardict.keys():
        msltd = get_mslt(rawfile)
        vardict["MSLT"].data = msltd
        print(f"Found variable {[vardict['MSLT'].varname, vardict['MSLT'].ramsname][rnameflag]}!")

    if "RTC" in vardict.keys():
        condmixd = get_condmix(vardict)
        vardict["RTC"].data = condmixd
        print(f"Found variable {[vardict['RTC'].varname, vardict['RTC'].ramsname][rnameflag]}!")
    
    if "RTI" in vardict.keys():
        icemixd = get_icemix(vardict)
        vardict["RTI"].data = icemixd
        print(f"Found variable {[vardict['RTI'].varname, vardict['RTI'].ramsname][rnameflag]}!")

    if "RTL" in vardict.keys():
        liqmixd = get_liqmix(vardict)
        vardict["RTL"].data = liqmixd
        print(f"Found variable {[vardict['RTL'].varname, vardict['RTL'].ramsname][rnameflag]}!")

    if "THV" in vardict.keys():
        thetavd = get_thetav(vardict)
        vardict["THV"].data = thetavd
        print(f"Found variable {[vardict['THV'].varname, vardict['THV'].ramsname][rnameflag]}!")

    if "THR" in vardict.keys():
        thetarhod = get_thetarho(vardict)
        vardict["THR"].data = thetarhod
        print(f"Found variable {[vardict['THR'].varname, vardict['THR'].ramsname][rnameflag]}!")

    if "P" in vardict.keys():
        pressured = get_pressure(vardict)
        vardict["P"].data = pressured
        print(f"Found variable {[vardict['P'].varname, vardict['P'].ramsname][rnameflag]}!")

    if "T" in vardict.keys():
        tempd = get_temp(vardict)
        vardict["T"].data = tempd
        print(f"Found variable {[vardict['T'].varname, vardict['T'].ramsname][rnameflag]}!")

    if "TV" in vardict.keys():
        tempvd = get_tempv(vardict)
        vardict["TV"].data = tempvd
        print(f"Found variable {[vardict['TV'].varname, vardict['TV'].ramsname][rnameflag]}!")

    if "TR" in vardict.keys():
        temprhod = get_temprho(vardict)
        vardict["TR"].data = temprhod
        print(f"Found variable {[vardict['TR'].varname, vardict['TR'].ramsname][rnameflag]}!")
    
    if "RHO" in vardict.keys():
        rhod = get_rho(vardict)
        vardict["RHO"].data = rhod
        print(f"Found variable {[vardict['RHO'].varname, vardict['RHO'].ramsname][rnameflag]}!")
    
    if "E" in vardict.keys():
        vaporpresd = get_vaporpres(vardict)
        vardict["E"].data = vaporpresd
        print(f"Found variable {[vardict['E'].varname, vardict['E'].ramsname][rnameflag]}!")

    if "PWAT" in vardict.keys():
        pwatd = get_pwat(rawfile, rawzsub)
        vardict["PWAT"].data = pwatd
        print(f"Found variable {[vardict['PWAT'].varname, vardict['PWAT'].ramsname][rnameflag]}!")
        
    if "TD" in vardict.keys():
        dewpointd = get_dewpoint(vardict)
        vardict["TD"].data = dewpointd
        print(f"Found variable {[vardict['TD'].varname, vardict['TD'].ramsname][rnameflag]}!")
    
    if "ZCT" in vardict.keys():
        cloudtopheightd = get_cloudtopheight(rawfile, rawzsub)
        vardict["ZCT"].data = cloudtopheightd
        print(f"Found variable {[vardict['ZCT'].varname, vardict['ZCT'].ramsname][rnameflag]}!")
    
    if "ZCB" in vardict.keys():
        cloudbaseheightd = get_cloudbaseheight(rawfile, rawzsub)
        vardict["ZCB"].data = cloudbaseheightd
        print(f"Found variable {[vardict['ZCB'].varname, vardict['ZCB'].ramsname][rnameflag]}!")

    if "PCT" in vardict.keys():
        cloudtoppresd = get_cloudtoppres(rawfile)
        vardict["PCT"].data = cloudtoppresd
        print(f"Found variable {[vardict['PCT'].varname, vardict['PCT'].ramsname][rnameflag]}!")

    if "PCB" in vardict.keys():
        cloudbasepresd = get_cloudbasepres(rawfile)
        vardict["PCB"].data = cloudbasepresd
        print(f"Found variable {[vardict['PCB'].varname, vardict['PCB'].ramsname][rnameflag]}!")
    
    #Vertical velocity advection
    if "WADVH" in vardict.keys() and ccoords is not None:
        wadvhd = get_wadv_h(vardict, gridprops)
        vardict["WADVH"].data = wadvhd
        print(f"Found variable {[vardict['WADVH'].varname, vardict['WADVH'].ramsname][rnameflag]}!")
    elif "WAVDH" in vardict.keys() and ccoords is None:
        print("Unfortunately, horizontal advection of w is only available for cartesian interpolation at the moment")
    
    
    if "WADVV" in vardict.keys() and ccoords is not None:
        wadvvd = get_wadv_v(vardict, ccoords)
        vardict["WADVV"].data = wadvvd
        print(f"Found variable {[vardict['WADVV'].varname, vardict['WADVV'].ramsname][rnameflag]}!")
    elif "WADVV" in vardict.keys() and ccoords is None:
        print("Unfortunately, vertical advection of w is only available for cartesian interpolation at the moment")

    #Vorticity
    if "VOX" in vardict.keys() and ccoords is not None:
        xvortd = get_xvort(vardict, gridprops, ccoords)
        vardict["VOX"].data = xvortd
        print(f"Found variable {[vardict['VOX'].varname, vardict['VOX'].ramsname][rnameflag]}!")
    elif "VOX" in vardict.keys() and ccoords is None:
        print("Unfortunately, x horizontal vorticity is only available for cartesian interpolation at the moment")

    if "VOY" in vardict.keys() and ccoords is not None:
        yvortd = get_yvort(vardict, gridprops, ccoords)
        vardict["VOY"].data = yvortd
        print(f"Found variable {[vardict['VOY'].varname, vardict['VOY'].ramsname][rnameflag]}!")
    elif "VOY" in vardict.keys() and ccoords is None:
        print("Unfortunately, y horizontal vorticity is only available for cartesian interpolation at the moment")

    if "VOZR" in vardict.keys() and ccoords is not None:
        zvortreld = get_zvortrel(vardict, gridprops)
        vardict["VOZR"].data = zvortreld
        print(f"Found variable {[vardict['VOZR'].varname, vardict['VOZR'].ramsname][rnameflag]}!")
    elif "VOZR" in vardict.keys() and ccoords is None:
        print("Unfortunately, relative vertical vorticity is only available for cartesian interpolation at the moment")
    
    if "VOZA" in vardict.keys() and ccoords is not None:
        zvortabsd = get_zvortabs(vardict, gridprops, ccoords)
        vardict["VOZA"].data = zvortabsd
        print(f"Found variable {[vardict['VOZA'].varname, vardict['VOZA'].ramsname][rnameflag]}!")
    elif "VOZA" in vardict.keys() and ccoords is None:
        print("Unfortunately, absolute vertical vorticity is only available for cartesian interpolation at the moment")
    
    #Supersaturation metrics
    if "SL" in vardict.keys():
        liqsatfracd = get_liqsatfrac(vardict)
        vardict["SL"].data = liqsatfracd
        print(f"Found variable {[vardict['SL'].varname, vardict['SL'].ramsname][rnameflag]}!")

    if "DVL" in vardict.keys():
        liqvapordeficitd = get_liqvapordeficit(vardict)
        vardict["DVL"].data = liqvapordeficitd
        print(f"Found variable {[vardict['DVL'].varname, vardict['DVL'].ramsname][rnameflag]}!")
    
    if "SI" in vardict.keys():
        icesatfracd = get_icesatfrac(vardict)
        vardict["SI"].data = icesatfracd
        print(f"Found variable {[vardict['SI'].varname, vardict['SI'].ramsname][rnameflag]}!")

    if "DVI" in vardict.keys():
        icevapordeficitd = get_icevapordeficit(vardict)
        vardict["DVI"].data = icevapordeficitd
        print(f"Found variable {[vardict['DVI'].varname, vardict['DVI'].ramsname][rnameflag]}!")

    if "RVSL" in vardict.keys():
        superliqmixd = get_superliqmix(vardict)
        vardict["RVSL"].data = superliqmixd
        print(f"Found variable {[vardict['RVSL'].varname, vardict['RVSL'].ramsname][rnameflag]}!")
    
    if "RVSI" in vardict.keys():
        supericemixd = get_supericemix(vardict)
        vardict["RVSI"].data = supericemixd
        print(f"Found variable {[vardict['RVSI'].varname, vardict['RVSI'].ramsname][rnameflag]}!")

    #Hydrometeor number concentrations
    if "NCP" in vardict.keys():
        cloudnumd = get_cloudnum(vardict)
        vardict["NCP"].data = cloudnumd
        print(f"Found variable {[vardict['NCP'].varname, vardict['NCP'].ramsname][rnameflag]}!")

    if "NDP" in vardict.keys():
        driznumd = get_drizzlenum(vardict)
        vardict["NDP"].data = driznumd
        print(f"Found variable {[vardict['NDP'].varname, vardict['NDP'].ramsname][rnameflag]}!")
    
    if "NRP" in vardict.keys():
        rainnumd = get_rainnum(vardict)
        vardict["NRP"].data = rainnumd
        print(f"Found variable {[vardict['NRP'].varname, vardict['NRP'].ramsname][rnameflag]}!")

    if "NPP" in vardict.keys():
        prisnumd = get_prisnum(vardict)
        vardict["NPP"].data = prisnumd
        print(f"Found variable {[vardict['NPP'].varname, vardict['NPP'].ramsname][rnameflag]}!")

    if "NSP" in vardict.keys():
        snownumd = get_snownum(vardict)
        vardict["NSP"].data = snownumd
        print(f"Found variable {[vardict['NSP'].varname, vardict['NSP'].ramsname][rnameflag]}!")

    if "NAP" in vardict.keys(): #I need a nap right now tbh
        aggnumd = get_aggregatenum(vardict)
        vardict["NAP"].data = aggnumd
        print(f"Found variable {[vardict['NAP'].varname, vardict['NAP'].ramsname][rnameflag]}!")

    if "NGP" in vardict.keys():
        graupelnumd = get_graupelnum(vardict)
        vardict["NGP"].data = graupelnumd
        print(f"Found variable {[vardict['NGP'].varname, vardict['NGP'].ramsname][rnameflag]}!")

    if "NHP" in vardict.keys():
        hailnumd = get_hailnum(vardict)
        vardict["NHP"].data = hailnumd
        print(f"Found variable {[vardict['NHP'].varname, vardict['NHP'].ramsname][rnameflag]}!")

    if "NIPP" in vardict.keys():
        platenumd = get_platenum(vardict)
        vardict["NIPP"].data = platenumd
        print(f"Found variable {[vardict['NIPP'].varname, vardict['NIPP'].ramsname][rnameflag]}!")
    
    if "NICP" in vardict.keys():
        columnnumd = get_columnnum(vardict)
        vardict["NICP"].data = columnnumd
        print(f"Found variable {[vardict['NICP'].varname, vardict['NICP'].ramsname][rnameflag]}!")

    if "NIDP" in vardict.keys():
        dendritenumd = get_dendritenum(vardict)
        vardict["NIDP"].data = dendritenumd
        print(f"Found variable {[vardict['NICP'].varname, vardict['NICP'].ramsname][rnameflag]}!")


    #Hydrometeor Water Contents
    if "CWC" in vardict.keys():
        cloudcond = get_cloudwatercontent(vardict)
        vardict["CWC"].data = cloudcond
        print(f"Found variable {[vardict['CWC'].varname, vardict['CWC'].ramsname][rnameflag]}!")
    
    if "DWC" in vardict.keys():
        drizcond = get_drizzlewatercontent(vardict)
        vardict["DWC"].data = drizcond
        print(f"Found variable {[vardict['DWC'].varname, vardict['DWC'].ramsname][rnameflag]}!")
    
    if "RWC" in vardict.keys():
        raincond = get_rainwatercontent(vardict)
        vardict["RWC"].data = raincond
        print(f"Found variable {[vardict['RWC'].varname, vardict['RWC'].ramsname][rnameflag]}!")
    
    if "PWC" in vardict.keys():
        priscond = get_priswatercontent(vardict)
        vardict["PWC"].data = priscond
        print(f"Found variable {[vardict['PWC'].varname, vardict['PWC'].ramsname][rnameflag]}!")
    
    if "SWC" in vardict.keys():
        snowcond = get_snowwatercontent(vardict)
        vardict["SWC"].data = snowcond
        print(f"Found variable {[vardict['SWC'].varname, vardict['SWC'].ramsname][rnameflag]}!")

    if "AWC" in vardict.keys():
        aggcond = get_aggregatewatercontent(vardict)
        vardict["AWC"].data = aggcond
        print(f"Found variable {[vardict['AWC'].varname, vardict['AWC'].ramsname][rnameflag]}!")
    
    if "GWC" in vardict.keys():
        graupelcond = get_graupelwatercontent(vardict)
        vardict["GWC"].data = graupelcond
        print(f"Found variable {[vardict['GWC'].varname, vardict['GWC'].ramsname][rnameflag]}!")

    if "HWC" in vardict.keys():
        hailcond = get_hailwatercontent(vardict)
        vardict["HWC"].data = hailcond
        print(f"Found variable {[vardict['HWC'].varname, vardict['HWC'].ramsname][rnameflag]}!")

    if "IPWC" in vardict.keys():
        platecond = get_platewatercontent(vardict)
        vardict["IPWC"].data = platecond
        print(f"Found variable {[vardict['IPWC'].varname, vardict['IPWC'].ramsname][rnameflag]}!")

    if "ICWC" in vardict.keys():
        columncond = get_columnwatercontent(vardict)
        vardict["ICWC"].data = columncond
        print(f"Found variable {[vardict['ICWC'].varname, vardict['ICWC'].ramsname][rnameflag]}!")

    if "IDWC" in vardict.keys():
        dendritecond = get_dendritewatercontent(vardict)
        vardict["IDWC"].data = dendritecond
        print(f"Found variable {[vardict['IDWC'].varname, vardict['IDWC'].ramsname][rnameflag]}!")
    

    #Vertically Integrated Hydrometeor Amounts
    if "VIC" in vardict.keys():
        cloudvertintd = get_cloudvertint(rawfile, rawzsub)
        vardict["VIC"].data = cloudvertintd
        print(f"Found variable {[vardict['VIC'].varname, vardict['VIC'].ramsname][rnameflag]}!")
    
    if "VID" in vardict.keys():
        drizvertintd = get_drizzlevertint(rawfile, rawzsub)
        vardict["VID"].data = drizvertintd
        print(f"Found variable {[vardict['VID'].varname, vardict['VID'].ramsname][rnameflag]}!")

    if "VIR" in vardict.keys():
        rainvertintd = get_rainvertint(rawfile, rawzsub)
        vardict["VIR"].data = rainvertintd
        print(f"Found variable {[vardict['VIR'].varname, vardict['VIR'].ramsname][rnameflag]}!")

    if "VIP" in vardict.keys():
        prisvertintd = get_prisvertint(rawfile, rawzsub)
        vardict["VIP"].data = prisvertintd
        print(f"Found variable {[vardict['VIP'].varname, vardict['VIP'].ramsname][rnameflag]}!")

    if "VIS" in vardict.keys():
        snowvertintd = get_snowvertint(rawfile, rawzsub)
        vardict["VIS"].data = snowvertintd
        print(f"Found variable {[vardict['VIS'].varname, vardict['VIS'].ramsname][rnameflag]}!")

    if "VIA" in vardict.keys():
        aggvertintd = get_aggregatevertint(rawfile, rawzsub)
        vardict["VIA"].data = aggvertintd
        print(f"Found variable {[vardict['VIA'].varname, vardict['VIA'].ramsname][rnameflag]}!")
    
    if "VIG" in vardict.keys():
        graupelvertintd = get_graupelvertint(rawfile, rawzsub)
        vardict["VIG"].data = graupelvertintd
        print(f"Found variable {[vardict['VIG'].varname, vardict['VIG'].ramsname][rnameflag]}!")

    if "VIH" in vardict.keys():
        hailvertintd = get_hailvertint(rawfile, rawzsub)
        vardict["VIH"].data = hailvertintd
        print(f"Found variable {[vardict['VIH'].varname, vardict['VIH'].ramsname][rnameflag]}!")

    if "VIPP" in vardict.keys():
        platevertintd = get_platevertint(rawfile, rawzsub)
        vardict["VIPP"].data = platevertintd
        print(f"Found variable {[vardict['VIIP'].varname, vardict['VIIP'].ramsname][rnameflag]}!")
    
    if "VICP" in vardict.keys():
        columnvertintd = get_columnvertint(rawfile, rawzsub)
        vardict["VICP"].data = columnvertintd
        print(f"Found variable {[vardict['VICP'].varname, vardict['VICP'].ramsname][rnameflag]}!")

    if "VIDP" in vardict.keys():
        dendritevertintd = get_dendritevertint(rawfile, rawzsub)
        vardict["VIDP"].data = dendritevertintd
        print(f"Found variable {[vardict['VIDP'].varname, vardict['VIDP'].ramsname][rnameflag]}!")

    if "VIL" in vardict.keys():
        vertintliqd = get_liqvertint(rawfile, rawzsub)
        vardict["VIL"].data = vertintliqd
        print(f"Found variable {[vardict['VIL'].varname, vardict['VIL'].ramsname][rnameflag]}!")

    if "VII" in vardict.keys():
        vertinticed = get_icevertint(rawfile, rawzsub)
        vardict["VII"].data = vertinticed
        print(f"Found variable {[vardict['VII'].varname, vardict['VII'].ramsname][rnameflag]}!")

    if "VIT" in vardict.keys():
        vertintcondd = get_condvertint(rawfile, rawzsub)
        vardict["VIT"].data = vertintcondd
        print(f"Found variable {[vardict['VIT'].varname, vardict['VIT'].ramsname][rnameflag]}!")

    #Hydrometeor Diameters
    if "DCP" in vardict.keys():
        clouddiamd = get_clouddiam(vardict, hydrodf)
        vardict["DCP"].data = clouddiamd
        print(f"Found variable {[vardict['DCP'].varname, vardict['DCP'].ramsname][rnameflag]}!")

    if "DDP" in vardict.keys():
        drizdiamd = get_drizzlediam(vardict, hydrodf)
        vardict["DDP"].data = drizdiamd
        print(f"Found variable {[vardict['DDP'].varname, vardict['DDP'].ramsname][rnameflag]}!")

    if "DRP" in vardict.keys():
        raindiamd = get_raindiam(vardict, hydrodf)
        vardict["DRP"].data = raindiamd
        print(f"Found variable {[vardict['DRP'].varname, vardict['DRP'].ramsname][rnameflag]}!")

    if "DPP" in vardict.keys():
        prisdiamd = get_prisdiam(vardict, hydrodf)
        vardict["DPP"].data = prisdiamd
        print(f"Found variable {[vardict['DPP'].varname, vardict['DPP'].ramsname][rnameflag]}!")

    if "DSP" in vardict.keys():
        snowdiamd = get_snowdiam(vardict, hydrodf)
        vardict["DSP"].data = snowdiamd
        print(f"Found variable {[vardict['DSP'].varname, vardict['DSP'].ramsname][rnameflag]}!")

    if "DAP" in vardict.keys():
        aggdiamd = get_aggregatediam(vardict, hydrodf)
        vardict["DAP"].data = aggdiamd
        print(f"Found variable {[vardict['DAP'].varname, vardict['DAP'].ramsname][rnameflag]}!")

    if "DGP" in vardict.keys():
        graupeldiamd = get_graupeldiam(vardict, hydrodf)
        vardict["DGP"].data = graupeldiamd
        print(f"Found variable {[vardict['DGP'].varname, vardict['DGP'].ramsname][rnameflag]}!")

    if "DHP" in vardict.keys():
        haildiamd = get_haildiam(vardict, hydrodf)
        vardict["DHP"].data = haildiamd
        print(f"Found variable {[vardict['DHP'].varname, vardict['DHP'].ramsname][rnameflag]}!")
    
    #Reference State and Momentum Budget Calculations. Currently, can only be done for cartesian interpolation, unfortunately. Be aware that these are *very slow* to calculate compared to all other 3D variables
    if "RHOP" in vardict.keys() and ccoords is not None:
        rhoprimed = get_rhoprime(vardict, window, kernname)
        vardict["RHOP"].data = rhoprimed
        print(f"Found variable {[vardict['RHOP'].varname, vardict['RHOP'].ramsname][rnameflag]}!")
    elif "RHOP" in vardict.keys() and ccoords is None:
        print("Unfortunately, perturbation density is only available for cartesian interpolation")
    
    if "BUOY_RHO" in vardict.keys() and ccoords is not None:
        rhobuoyd = get_rhobuoy(vardict, window, kernname)
        vardict["BUOY_RHO"].data = rhobuoyd
        print(f"Found variable {[vardict['BUOY_RHO'].varname, vardict['BUOY_RHO'].ramsname][rnameflag]}!")
    elif "BUOY_RHO" in vardict.keys() and ccoords is None: 
        print("Unfortunately, buoyancy calculations are only available for cartesian interpolation")
    
    if "BUOY_THETA" in vardict.keys() and ccoords is not None:
        thetabuoyd = get_thetabuoy(vardict, window, kernname)
        vardict["BUOY_THETA"].data = thetabuoyd
        print(f"Found variable {[vardict['BUOY_THETA'].varname, vardict['BUOY_THETA'].ramsname][rnameflag]}!")
    elif "BUOY_THETA" in vardict.keys() and ccoords is None:
        print("Unfortunately, buoyancy calculations are only available for cartesian interpolation")
    
    if "BUOY_TEMP" in vardict.keys() and ccoords is not None:
        thermbuoyd = get_thermalbuoy(vardict, window, kernname)
        vardict["BUOY_TEMP"].data = thermbuoyd
        print(f"Found variable {[vardict['BUOY_TEMP'].varname, vardict['BUOY_TEMP'].ramsname][rnameflag]}!")
    elif "BUOY_TEMP" in vardict.keys() and ccoords is None:
        print("Unfortunately, buoyancy calculations are only available for cartesian interpolation")

    if "BUOY_VAP" in vardict.keys() and ccoords is not None:
        vapbuoyd = get_vaporbuoy(vardict, window, kernname)
        vardict["BUOY_VAP"].data = vapbuoyd
        print(f"Found variable {[vardict['BUOY_VAP'].varname, vardict['BUOY_VAP'].ramsname][rnameflag]}!")
    elif "BUOY_VAP" in vardict.keys() and ccoords is None:
        print("Unfortunately, buoyancy calculations are only available for cartesian interpolation")

    if "BUOY_COND" in vardict.keys() and ccoords is not None:
        condbuoyd = get_condbuoy(vardict)
        vardict["BUOY_COND"].data = condbuoyd
        print(f"Found variable {[vardict['BUOY_COND'].varname, vardict['BUOY_COND'].ramsname][rnameflag]}!")
    elif "BUOY_COND" in vardict.keys() and ccoords is None:
        print("Unfortunately, buoyancy calculations are only available for cartesian interpolation")
    
    if "BUOY_PPRIME" in vardict.keys() and ccoords is not None:
        ppbuoyd = get_pprimebuoy(vardict, window, kernname)
        vardict["BUOY_PPRIME"].data = ppbuoyd
        print(f"Found variable {[vardict['BUOY_PPRIME'].varname, vardict['BUOY_PPRIME'].ramsname][rnameflag]}!")
    elif "BUOY_PPRIME" in vardict.keys() and ccoords is None:
        print("Unfortunately, buoyancy calculations are only available for cartesian interpolation") 
    
    if "VPPGF" in vardict.keys() and ccoords is not None:
        vppgfd = get_vppgf(vardict, window, kernname, ccoords)
        vardict["VPPGF"].data = vppgfd
        print(f"Found variable {[vardict['VPPGF'].varname, vardict['VPPGF'].ramsname][rnameflag]}!")
    elif "VPPGF" in vardict.keys() and ccoords is None:
        print("Unfortunately, VPPGF is only available for cartesian interpolation") 
    return vardict