#!/usr/bin/env python
# coding: utf-8

# In[2]:

import pandas as pd
from varinit import outvar
from interplib import get_zact
from derivcalcs import *
from dvardict import derivvarinit
import warnings

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
CloudWaterContent: CWC --> Volumetric water content of cloud droplets, in kg/m^3
DrizzleWaterContent: DWC --> Volumetric water content of drizzle, in kg/m^3
RainWaterContent: RWC --> Volumetric water content of rain, in kg/m^3
PrisWaterContent: PWC --> Volumetric water content of pristine ice, in kg/m^3
SnowWaterContent: SWC --> Volumetric water content of snow, in kg/m^3
AggWaterContent: AWC --> Volumetric water content of aggregates, in kg/m^3
GraupelWaterContent: GWC --> Volumetric water content of graupel, in kg/m^3
HailWaterContent: HWC, --> Volumetric water content of hail, in kg/m^3
PlateWaterContent: IPWC --> Volumetric water content of plate ice, in kg/m^3
ColumnWaterContent: ICWC --> Volumetric water content of column ice, in kg/m^3
DendriteWaterContent: IDWC --> Volumetric water content of dendrites, in kg/m^3
VertIntCloud: VIC --> Mass per unit area of all cloud droplets in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
VertIntDrizzle: VID --> Mass per unit area of all drizzle in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
VertIntRain: VIR --> Mass per unit area of all rain in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
VertIntPris: VIP --> Mass per unit area of all pristine ice in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
VertIntSnow: VIS --> Mass per unit area of all snow in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
VertIntAgg: VIA --> Mass per unit area of all aggregates in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
VertIntGraupel: VIG --> Mass per unit area of all graupel in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
VertIntHail: VIH --> Mass per unit area of all hail in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
VertIntPlate: VIIP --> Mass per unit area of all plate ice in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
VertIntColumn: VIIC --> Mass per unit area of all column ice in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
VertIntDendrite: VIID --> Mass per unit area of all dendrites in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
VertIntLiq: VIL --> Mass per unit area of total liquid condensate in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
VertIntIce: VII --> Mass per unit area of total ice condensate in the column. Natural units are kg/m^2, but I've converted that to a liquid water equivalent depth in mm (like PWAT)
CloudDiam: DCP --> Diameter of cloud droplets, in mm
DrizzleDiam: DDP --> Diameter of drizzle, in mm
RainDiam: DRP --> Diameter of rain, in mm
PrisDiam: DPP --> Diameter of pristine ice, in mm
SnowDiam: DSP --> Diameter of snow, in mm
AggDiam: DAP --> Diameter of aggregates, in mm
GraupelDiam: DGP --> Diameter of graupel, in mm
HailDiam: DHP --> Diameter of hail, in mm
RhoPrime: RHOP --> Perturbation density (current grid cell density subtracted from horizontal mean density calculated based on convolution kernel and window size), in kg/m^3
RhoBuoy: BUOY_RHO --> Total buoyancy (buoyancy vertical acceleration) calculated from perturbation density, in m/s^2. Has a larger residual than ThetaRhoBuoy
ThetaRhoBuoy: BUOY_THETA --> Total buoyancy (buoyancy vertical acceleration) calculated from perturbation density potential temperature, in m/s^2. Has a smaller residual than RhoBuoy
ThermalBuoy: BUOY_TEMP --> Vertical acceleration from thermal buoyancy, in m/s^2.
VapBuoy: BUOY_VAP --> Vertical acceleration from vapor buouancy, in m/s^2.
CondBuoy: BUOY_COND --> Vertical acceleration from condensate loading, in m/s^2.
PPrimeBuoy: BUOY_PPRIME --> Vertical acceleration from horizontal pressure gradients, in m/s^2.
VPPGF: VPPGF --> Vertical acceleration from vertical pressure gradients (Vertical Perturbation Pressure Gradient Force), in m/s^2.
'''

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
    with warnings.catch_warnings(): #The reason we're silenceing warnings here is because numpy get's annoyed if a slice is all NaNs. Since we're looking at cloud tops and bottoms, any column without any condensate is going to be all NaNs!
        warnings.filterwarnings("ignore")
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
    with warnings.catch_warnings(): #The reason we're silencing warnings here is because the diameter calculation involves division by particle number mixing ratio. Obviously, if there are no particles of that species, this will be a divison by zero and numpy will throw a fit. Since we don't care about the diameter of particles in grid cells without any particles, we can silence the warning.
        warnings.filterwarnings("ignore")
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