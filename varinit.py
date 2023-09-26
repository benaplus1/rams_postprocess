import xarray as xr
import numpy as np
from time import perf_counter
from numpy import float32, float64
# from functools import partial
# from concurrent.futures import ProcessPoolExecutor

class outvar:
    def __init__(self, varname = "", longname = "", stdname = None, vector = False, offline = False, vartype = "", nptype = float32, ramsname = "", invar = None, unitfactor = "", decnum = 3, units = "", data = None):
        #varname (string) is the name of the variable as referred to in the xarray dataset
        #longname (string) is the "long name" of the variable as shown in the xarray attributes
        #stdname (string) is the "standard name" according to the cf naming conventions (if applicable), to be used in the xarray attributes
        #vector is whether this variable is a vector (u,v,w) that needs a special interpolation routine
        #offline (bool) is whether this variable is simply a renamed and truncated version of a raw RAMS variable (ex. Theta, RainPrecipRate) or an offline-calculated variable (like Pressure or CloudTopHeight)
        #vartype (string) is whether this variable is 3d (z,y,x), 2d (y,x), 2d-patch (patch, y, x), 2d-snow (patch, snowlev, y, x), or 2d-soil (patch, soillev, y, x)
        #ntype (numpy dtype) is the dtype to use (float32 or float64) for the variable
        #ramsname (string) is the name of the raw RAMS variable read in before truncation and renaming. Derived variables have a "ramsname" which is in a similar style to the native RAMS variables
        #invar (list of strings) is the other variables (if any) that need to be passed in to calculate this variable. Only applies for offline-calculated variables
        #decnum (int) is the number of decimal points to use for storing this variable. This will only come into play after all unit conversions
        #unitfactor (float or string) is the factor by which this variable needs to be modified by its original units (for example, 2.5*10^6 for SFLUX_R (kg/m^2*s) to LatHeatFlux (W/m^2) or 1000 for RCP (kg/kg) to CloudMassMix (g/kg), 1000/afiletime for AGGREGATET to Aggregation (change in mixing ratio due to aggregation in the previous timestep)). These unit conversions are only done after all offline variables have been calculated
        #units (string) are the units to be put in the xarray attributes

        self.varname = varname
        self.longname = longname
        self.stdname = stdname
        self.vector = vector
        self.offline = offline
        self.vartype = vartype
        self.nptype = nptype
        self.ramsname = ramsname
        self.invar = invar
        self.unitfactor = unitfactor
        self.decnum = decnum
        self.units = units
        self.data = data #Assigned later in postprocess.py

    def __str__(self):
        return(f"Variable: {self.longname}. CF Standard Name: {self.stdname}. RAMS-Format Name: {self.ramsname}; Verbose Name: {self.varname}; Offline Variable: {self.offline}; Variable Dimensions: {self.vartype}; Vector Variable: {self.vector}; Other Variables Used to Calculate it: {self.invar}; Factor by Which This Variable is Multiplied at the End to Get Final Units: {self.unitfactor}; Number of Decimals Used to Store Variable: {self.decnum}; Units of Variable: {self.units};")
def fillvarlist():
    #Fullvarlist initializes a list of all RAMS variables which can be output (except KPP variables, which aren't implemented yet). The function below, get_vardict, actually figures out which variables are available and the user wants to include. The post-processing routine assumes that you want the same variables for all times in the post-processing. If you want to analyze new variables at certain times, restart the program with a new starttime and a new variable list file.
    varlist = []
    #3D Vector variables
    varlist.append(outvar(varname = "u", longname = "Eastward Wind Velocity", stdname = "eastward_wind", vector = True, vartype = "3d", nptype = float32, ramsname = "UC", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))
    varlist.append(outvar(varname = "v", longname = "Northward Wind Velocity", stdname = "northward_wind", vector = True, vartype = "3d", nptype = float32, ramsname = "VC", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))
    varlist.append(outvar(varname = "w", longname = "Vertical Velocity", stdname = "upward_air_velocity", vector = True, vartype = "3d", nptype = float32, ramsname = "WC", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))
    varlist.append(outvar(varname = "PastU", longname = "Past Timestep Eastward Wind Velocity", stdname = None, vector = True, vartype = "3d", nptype = float32, ramsname = "UP", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))
    varlist.append(outvar(varname = "PastV", longname = "Past Timestep Northward Wind Velocity", stdname = None, vector = True, vartype = "3d", nptype = float32, ramsname = "VP", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))
    varlist.append(outvar(varname = "PastW", longname = "Past Timestep Vertical Velocity", stdname = None, vector = True, vartype = "3d", nptype = float32, ramsname = "WP", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))

#3D Meteorology Scalar Variables
    varlist.append(outvar(varname = "VaporMix", longname = "Vapor Mixing Ratio", stdname = "humidity_mixing_ratio", vartype = "3d", nptype = float32, ramsname = "RV", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "WaterMix", longname = "Water Mixing Ratio", stdname = "water_mixing_ratio", vartype = "3d", nptype = float32, ramsname = "RTP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "Theta", longname = "Potential Temperature", stdname = "air_potential_temperature", vartype = "3d", nptype = float32, ramsname = "THETA", invar = None, unitfactor = 1, decnum = 3, units = "K"))
    varlist.append(outvar(varname = "Exner", longname = "Exner Function * Cp", stdname = None, vartype = "3d", nptype = float32, ramsname = "PI", invar = None, unitfactor = 1, decnum = 3, units = 1))
    varlist.append(outvar(varname = "ReferenceDensity", longname = "Reference Air Density", stdname = None, vartype = "3d", nptype = float32, ramsname = "DN0", invar = None, unitfactor = 1, decnum = 3, units = "kg m**-3"))
    varlist.append(outvar(varname = "ReferenceExner", longname = "Reference Exner Function * Cp", stdname = None, vartype = "3d", nptype = float32, ramsname = "PI0", invar = None, unitfactor = 1, decnum = 3, units = 1))
    varlist.append(outvar(varname = "ReferenceTheta", longname = "Reference Potential Temperature", stdname = None, vartype = "3d", nptype = float32, ramsname = "TH0", invar = None, unitfactor = 1, decnum = 3, units = "K"))
    varlist.append(outvar(varname = "PastExnerPrime", longname = "Past Perturbation Exner Function * Cp", stdname = None, vartype = "3d", ramsname = "PP", invar = None, unitfactor = 1, decnum = 3, units = 1))
    varlist.append(outvar(varname = "ExnerPrime", longname = "Current Perturbation Exner Function * Cp", stdname = None, vartype = "3d", ramsname = "PC", invar = None, unitfactor = 1, decnum = 3, units = 1))
    varlist.append(outvar(varname = "Thetail", longname = "Ice-Liquid Potential Temperature", stdname = None, vartype = "3d", ramsname = "THP", invar = None, unitfactor = 1, decnum = 3, units = "K"))

#3D Turbulence
    varlist.append(outvar(varname = "TKE", longname = "Turbulent Kinetic Energy", stdname = "specific_turbulent_kinetic_energy_of_air", vartype = "3d", ramsname = "TKEP", invar = None, unitfactor = 1, decnum = 3, units = "J kg**-1"))
    varlist.append(outvar(varname = "HorizHeatDiffusivityUnWeighted", longname = "Horizontal Eddy Diffusivity Coefficient for Heat (No Density Weighting)", stdname = None, vartype = "3d", ramsname = "HKH", invar = None, unitfactor = 1, decnum = 3, units = "m**2 s**-1"))
    varlist.append(outvar(varname = "VertHeatDiffusivityUnWeighted", longname = "Vertical Eddy Diffusivity Coefficient for Heat (No Density Weighting)", stdname = None, vartype = "3d", ramsname = "VKH", invar = None, unitfactor = 1, decnum = 3, units = "m**2 s**-1"))
    varlist.append(outvar(varname = "VertHeatDiffusivityWeighted", longname = "Vertical Eddy Diffusivity Coefficient for Heat", stdname = "atmosphere_heat_diffusivity", vartype = "3d", ramsname = "RVKH", invar = None, unitfactor = 1, decnum = 3, units = "m**2 s**-1"))
    varlist.append(outvar(varname = "HorizMomentumDiffusivityWeighted", longname = "Horizontal Eddy Diffusivity Coefficient for Momentum", stdname = None, vartype = "3d", ramsname = "RHKM", invar = None, unitfactor = 1, decnum = 3, units = "m**2 s**-1"))
    varlist.append(outvar(varname = "VertMomentumDiffusivityWeighted", longname = "Vertical Eddy Diffusivity Coefficient for Momentum", stdname = "atmosphere_momentum_diffusivity", vartype = "3d", ramsname = "RVKM", invar = None, unitfactor = 1, decnum = 3, units = "m**2 s**-1"))

#3d Convective Parameterization Variables
    varlist.append(outvar(varname = "ConvectiveHeating", longname = "Rate of Change in Air Temperature from Convective Parameterization", stdname = "tendency_of_air_temperature_due_to_convection", vartype = "3d", ramsname = "THSRC", invar = None, unitfactor = 1, decnum = 4, units = "K s**-1"))
    varlist.append(outvar(varname = "ConvectiveMoistening", longname = "Rate of Change in Vapor Mixing Ratio from Convective Parameterization", stdname = None, vartype = "3d", ramsname = "RTSRC", unitfactor = "gperkg", decnum = 4, units = "g kg**-1 s**-1"))
    varlist.append(outvar(varname = "ConvectiveCloudMixChange", longname = "Rate of Change in Cloud Droplet Mixing Ratio from Convective Parameterization", stdname = None, vartype = "3d", ramsname = "RCSRC", unitfactor = 10**6, decnum = 4, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "ConvectiveRainMixChange", longname = "Rate of Change in Rain Mixing Ratio from Convective Parameterization", stdname = None, vartype = "3d", ramsname = "RRSRC", unitfactor = 10**6, decnum = 4, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "ConvectivePrisMixChange", longname = "Rate of Change in Pristine Ice Mixing Ratio from Convective Parameterization", stdname = None, vartype = "3d", ramsname = "RPSRC", unitfactor = 10**6, decnum = 4, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "ConvectiveSnowMixChange", longname = "Rate of Change in Snow Mixing Ratio from Convective Parameterization", stdname = None, vartype = "3d", ramsname = "RSSRC", unitfactor = 10**6, decnum = 4, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "ConvectiveVertVelocity", longname = "Running Mean Average of Vertical Velocity", stdname = None, vartype = "3d", ramsname = "W0AVG", unitfactor = 1, decnum = 3, units = "m s**-1"))
    varlist.append(outvar(varname = "ConvectiveVertVelocityContra", longname = "Running Mean Average of Horizontal Components of Contravariant Vertical Velocity", stdname = None, vartype = "3d", ramsname = "W0AVGLT", unitfactor = 1, decnum = 3, units = "m s**-1"))

#3D Radiation Variables
    varlist.append(outvar(varname = "VisLength", longname = "Length of Visibility", stdname = "visibility_in_air", vartype = "3d", nptype = float32, ramsname = "BEXT", invar = None, unitfactor = 1000, decnum = 3, units = "m"))
    varlist.append(outvar(varname = "RadiativeHeating", longname = "Radiative Heating Rate in Air", stdname = "tendency_of_air_temperature_due_to_radiative_heating", vartype = "3d", nptype = float32, ramsname = "FTHRD", invar = None, unitfactor = 1, decnum = 3, units = "K s**-1"))
    varlist.append(outvar(varname = "LongwaveDown", longname = "Downwelling Longwave Radiative Flux", stdname = "downwelling_longwave_flux_in_air", vartype = "3d", nptype = float32, ramsname = "LWDN", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "LongwaveUp", longname = "Upwelling Longwave Radiative Flux", stdname = "upwelling_longwave_flux_in_air", vartype = "3d", nptype = float32, ramsname = "LWUP", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "ShortwaveDown", longname = "Downwelling Shortwave Radiative Flux", stdname = "downwelling_shortwave_flux_in_air", vartype = "3d", nptype = float32, ramsname = "SWDN", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "ShortwaveUp", longname = "Upwelling Shortave Radiative Flux", stdname = "upwelling_shortwave_flux_in_air", vartype = "3d", nptype = float32, ramsname = "SWUP", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))

#Hydrometeor Mass/Number Mixing Ratios
    varlist.append(outvar(varname = "CloudNumberMix", longname = "Cloud Droplet Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CCP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    varlist.append(outvar(varname = "DrizzleNumberMix", longname = "Drizzle Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CDP", invar = None, unitfactor = 1, decnum = 0, units = "# kg**-1"))
    varlist.append(outvar(varname = "RainNumberMix", longname = "Rain Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CRP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    varlist.append(outvar(varname = "SnowNumberMix", longname = "Snow Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CSP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    varlist.append(outvar(varname = "PrisNumberMix", longname = "Pristine Ice Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CPP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    varlist.append(outvar(varname = "AggNumberMix", longname = "Aggregate Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CAP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    varlist.append(outvar(varname = "GraupelNumberMix", longname = "Graupel Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CGP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    varlist.append(outvar(varname = "HailNumberMix", longname = "Hail Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CHP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    varlist.append(outvar(varname = "PlateNumberMix", longname = "Plate Ice Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CIPP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    varlist.append(outvar(varname = "ColumnIceNumberMix", longname = "Column Ice Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CICP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    varlist.append(outvar(varname = "DendriteNumberMix", longname = "Dendritic Snowflake Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CIDP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    varlist.append(outvar(varname = "CloudMassMix", longname = "Cloud Droplet Mixing Ratio", stdname = "cloud_liquid_water_mixing_ratio", vartype = "3d", nptype = float32, ramsname = "RCP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "DrizzleMassMix", longname = "Drizzle Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RDP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "RainMassMix", longname = "Rain Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RRP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "SnowMassMix", longname = "Snow Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RSP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "PrisMassMix", longname = "Pristine Ice Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RPP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "AggMassMix", longname = "Aggregate Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RAP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "GraupelMassMix", longname = "Graupel Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RGP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "HailMassMix", longname = "Hail Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RHP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "PlateMassMix", longname = "Plate Ice Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RIPP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "ColumnIceMassMix", longname = "Column Ice Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RICP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "DendriteMassMix", longname = "Dendritic Snowflake Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RIDP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))

#Hydrometeor Energies
    varlist.append(outvar(varname = "RainEnergy", longname = "Internal Energy of Rain (J/kg Hydrometeor)", stdname = None, vartype = "3d", ramsname = "Q2", invar = None, unitfactor = 1, decnum = 3, units = "J kg**-1"))
    varlist.append(outvar(varname = "GraupelEnergy", longname = "Internal Energy of Graupel (J/kg Hydrometeor)", stdname = None, vartype = "3d", ramsname = "Q6", invar = None, unitfactor = 1, decnum = 3, units = "J kg**-1"))
    varlist.append(outvar(varname = "HailEnergy", longname = "Internal Energy of Hail (J/kg Hydrometeor)", stdname = None, vartype = "3d", ramsname = "Q7", invar = None, unitfactor = 1, decnum = 3, units = "J kg**-1"))

#Hydrometeor growth budgets
    varlist.append(outvar(varname = "Aggregation", longname = "Rate of Aggregation Between Analysis File Times", stdname = None, vartype = "3d", nptype = float32, ramsname = "AGGREGATET", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "AggSnowPris", longname = "Rate of Aggregation from Collisions Between Pristine Ice and Snow Between Analysis File Times", stdname = None, vartype = "3d", nptype = float32, ramsname = "AGGRPRISSNOWT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "AggSelfSnow", longname = "Rate of Aggregation from Collisions Between Aggregates and Snow Between Analysis File Times", stdname = None, vartype = "3d", nptype = float32, ramsname = "AGGRSELFSNOWT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "AggSelfPris", longname = "Rate of Aggregation from Collisions Between Aggregates and Pristine Ice Between Analysis File Times", stdname = None, vartype = "3d", nptype = float32, ramsname = "AGGRSELFPRIST", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "LatentHeatFreezing", longname = "Rate of Temperature Change from Freezing/Melting", stdname = None, vartype = "3d", ramsname = "LATHEATFRZ", invar = None, unitfactor = "pertimestep", decnum = 6, units = "K s**-1"))
    varlist.append(outvar(varname = "LatentHeatFreezingTotal", longname = "Total Temperature Change from Freezing/Melting Between Analysis File Times", stdname = None, vartype = "3d", ramsname = "LATHEATFRZT", invar = None, unitfactor = 1, decnum = 3, units = "K"))
    varlist.append(outvar(varname = "LatentHeatCondensation", longname = "Rate of Temperature Change from Condensation/Evaporation and Deposition/Sublimation", stdname = None, vartype = "3d", ramsname = "LATHEATVAP", invar = None, unitfactor = "pertimestep", decnum = 6, units = "K s**-1"))
    varlist.append(outvar(varname = "LatentHeatCondensationTotal", longname = "Total Temperature Change from Condensation/Evaporation and Deposition/Sublimation Between Analysis File Times", stdname = None, vartype = "3d", ramsname = "LATHEATVAPT", invar = None, unitfactor = 1, decnum = 3, units = "K"))
    varlist.append(outvar(varname = "CollectedCloudWater", longname = "Rate of Collection of Cloud Droplets by Rain", stdname = None, vartype = "3d", ramsname = "CLD2RAINT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "IceCollectionMelting", longname = "Rate of Ice Melting from Collection of Rain", stdname = None, vartype = "3d", ramsname = "ICE2RAINT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "PrisRainFreezing", longname = "Rate of Ice Growth from Collisions of Rain and Pristine Ice", stdname = None, vartype = "3d", ramsname = "RAIN2PRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "SnowRainRiming", longname = "Rate of Snow Growth from Riming of Rain", stdname = None, vartype = "3d", ramsname = "RAIN2SNT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    varlist.append(outvar(varname = "AggRainRiming", longname = "Rate of Aggregate Growth from Riming of Rain", stdname = None, vartype = "3d", ramsname = "RAIN2AGT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "GraupelRainRiming", longname = "Rate of Graupel Growth from Riming of Rain", stdname = None, vartype = "3d", ramsname = "RAIN2GRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "HailRainRiming", longname = "Rate of Hail Growth from Riming of Rain", stdname = None, vartype = "3d", ramsname = "RAIN2HAT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "SnowCloudRiming", longname = "Rate of Snow Growth from Riming of Cloud Droplets", stdname = None, vartype = "3d", ramsname = "RIMECLDSNOWT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "AggCloudRiming", longname = "Rate of Aggregate Growth from Riming of Cloud Droplets", stdname = None, vartype = "3d", ramsname = "RIMECLDAGGRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "GraupelCloudRiming", longname = "Rate of Graupel Growth from Riming of Cloud Droplets", stdname = None, vartype = "3d", ramsname = "RIMECLDGRAUT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "HailCloudRiming", longname = "Rate of Hail Growth from Riming of Cloud Droplets", stdname = None, vartype = "3d", ramsname = "RIMECLDHAILT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "TotalCloudRiming", longname = "Rate of Total Ice Growth from Riming of Cloud Droplets", stdname = None, vartype = "3d", ramsname = "RIMECLDT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "TotalRainRiming", longname = "Rate of Total Ice Growth from Riming of Rain", stdname = None, vartype = "3d", ramsname = "RAIN2ICET", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "CloudVapGrowth", longname = "Rate of Cloud Droplet Growth from Condensation", stdname = None, vartype = "3d", ramsname = "VAPCLDT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "DrizzleVapGrowth", longname = "Rate of Drizzle Growth from Condensation", stdname = None, vartype = "3d", ramsname = "VAPDRIZT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "RainVapGrowth", longname = "Rate of Rain Growth from Condensation", stdname = None, vartype = "3d", ramsname = "VAPRAINT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "TotalLiquidVapGrowth", longname = "Rate of Total Liquid Growth from Condensation", stdname = None, vartype = "3d", ramsname = "VAPLIQT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "PrisVapGrowth", longname = "Rate of Pristine Ice Growth from Vapor Deposition", stdname = None, vartype = "3d", ramsname = "VAPPRIST", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    varlist.append(outvar(varname = "SnowVapGrowth", longname = "Rate of Snow Growth from Vapor Deposition", stdname = None, vartype = "3d", ramsname = "VAPSNOWT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    varlist.append(outvar(varname = "AggVapGrowth", longname = "Rate of Aggregate Growth from Vapor Deposition", stdname = None, vartype = "3d", ramsname = "VAPAGGRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    varlist.append(outvar(varname = "GraupelVapGrowth", longname = "Rate of Graupel Growth from Vapor Deposition", stdname = None, vartype = "3d", ramsname = "VAPGRAUT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    varlist.append(outvar(varname = "HailVapGrowth", longname = "Rate of Hail Growth from Vapor Deposition", stdname = None, vartype = "3d", ramsname = "VAPHAILT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    varlist.append(outvar(varname = "TotalIceVapGrowth", longname = "Rate of Total Ice Growth from Vapor Deposition", stdname = None, vartype = "3d", ramsname = "VAPICET", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    varlist.append(outvar(varname = "CloudEvap", longname = "Rate of Cloud Droplet Mass Loss from Evaporation", stdname = None, vartype = "3d", ramsname = "EVAPCLDT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "DrizzleEvap", longname = "Rate of Drizzle Mass Loss from Evaporation", stdname = None, vartype = "3d", ramsname = "EVAPCLDT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "RainEvap", longname = "Rate of Rain Mass Loss from Evaporation", stdname = None, vartype = "3d", ramsname = "EVAPRAINT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "TotalLiquidEvap", longname = "Rate of Total Liquid Water Mass Loss from Evaporation", stdname = None, vartype = "3d", ramsname = "EVAPLIQT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "PristineEvap", longname = "Rate of Pristine Ice Mass Loss from Sublimation", stdname = None, vartype = "3d", ramsname = "EVAPPRIST", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "SnowEvap", longname = "Rate of Snow Mass Loss from Sublimation", stdname = None, vartype = "3d", ramsname = "EVAPPRIST", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "AggEvap", longname = "Rate of Aggregate Mass Loss from Sublimation", stdname = None, vartype = "3d", ramsname = "EVAPAGGRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "GraupelEvap", longname = "Rate of Graupel Mass Loss from Sublimation", stdname = None, vartype = "3d", ramsname = "EVAPGRAUT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "HailEvap", longname = "Rate of Hail Mass Loss from Sublimation", stdname = None, vartype = "3d", ramsname = "EVAPHAILT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "PrisMelting", longname = "Rate of Pristine Ice Melting", stdname = None, vartype = "3d", ramsname = "MELTPRIST", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "SnowMelting", longname = "Rate of Snow Melting", stdname = None, vartype = "3d", ramsname = "MELTSNOWT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "AggMelting", longname = "Rate of Aggregate Melting", stdname = None, vartype = "3d", ramsname = "MELTAGGRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "GraupelMelting", longname = "Rate of Graupel Melting", stdname = None, vartype = "3d", ramsname = "MELTGRAUT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "HailMelting", longname = "Rate of Hail Melting", stdname = None, vartype = "3d", ramsname = "MELTHAILT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "TotalIceMelting", longname = "Rate of Total Ice Melting", stdname = None, vartype = "3d", ramsname = "MELTICET", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))

#Freezing Processes
    varlist.append(outvar(varname = "ContactFreezing", longname = "Rate of Contact Freezing", stdname = None, vartype = "3d", ramsname = "INUCCONTRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "HomogeneousFreezing", longname = "Rate of Homogeneous Freezing", stdname = None, vartype = "3d", ramsname = "INUCHOMRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "ImmersionFreezing", longname = "Rate of Condensation/Immersion Freezing", stdname = None, vartype = "3d", ramsname = "INUCIFNRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "HazeFreezing", longname = "Rate of Haze Droplet Freezing", stdname = None, vartype = "3d", ramsname = "INUCHAZRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))

#Dust Nucleation
    varlist.append(outvar(varname = "Dust1CloudNucleation", longname = "Rate of Nucleation of Cloud Droplets from Dust Mode 1", stdname = None, vartype = "3d", ramsname = "DUST1CLDRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "Dust2CloudNucleation", longname = "Rate of Nucleation of Cloud Droplets from Dust Mode 2", stdname = None, vartype = "3d", ramsname = "DUST2CLDRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "Dust1DrizzleNucleation", longname = "Rate of Nucleation of Drizzle from Dust Mode 1", stdname = None, vartype = "3d", ramsname = "DUST1DRZRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    varlist.append(outvar(varname = "Dust2DrizzleNucleation", longname = "Rate of Nucleation of Drizzle from Dust Mode 2", stdname = None, vartype = "3d", ramsname = "DUST2CLDRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))

#Aerosol Concentrations and Processes
    varlist.append(outvar(varname = "CCNNumberMix", longname = "Cloud Condensation Nuclei Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CCCNP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "INPNumberMix", longname = "Ice Nucleating Particle Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CIFNP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "GCCNNumberMix", longname = "Giant CCN Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "GCCNP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "CN1NumberMix", longname = "CCN Mode 1 Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CN1NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "CN2NumberMix", longname = "CCN Mode 2 Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CN2NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "SubDustNumberMix", longname = "Sub-Micron Dust Number Concentration", stdname = None, vartype = "3d", ramsname = "MD1NP", invar = None, unitfactor = 10**-6, decnum = 5, units = "# mg**-1"))
    varlist.append(outvar(varname = "SuperDustNumberMix", longname = "Super-Micron Dust Number Concentration", stdname = None, vartype = "3d", ramsname = "MD2NP", invar = None, unitfactor = 10**-6, decnum = 5, units = "# mg**-1"))
    varlist.append(outvar(varname = "ABC1NumberMix", longname = "Absorbing Carbon Type 1 Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "ABC1NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "ABC2NumberMix", longname = "Absorbing Carbon Type 2 Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "ABC2NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "SaltFilmNumberMix", longname = "Sea Salt Film Drop Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "SALT_FILM_NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "SaltJetNumberMix", longname = "Sea Salt Jet Drop Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "SALT_JET_NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "SaltSpumeNumberMix", longname = "Sea Salt Spume Drop Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "SALT_SPUM_NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "SubAeroRegenNumberMix", longname = "Sub-Micron Regenerated Aerosol Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "REGEN_AERO1_NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "SuperAeroRegenNumberMix", longname = "Super-Micron Regenerated Aerosol Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "REGEN_AERO2_NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "CCNMassMix", longname = "Cloud Condensation Nuclei Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CCCMP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "GCCNMassMix", longname = "Giant CCN Mixing Ratio", stdname = None, vartype = "3d", ramsname = "GCCMP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "CN1MassMix", longname = "CCN Mode 1 Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CN1MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "CN2MassMix", longname = "CCN Mode 2 Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CN2MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "SubDustMassMix", longname = "Sub-Micron Dust Mixing Ratio", stdname = None, vartype = "3d", ramsname = "MD1MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "SuperDustMassMix", longname = "Super-Micron Dust Mixing Ratio", stdname = None, vartype = "3d", ramsname = "MD2MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "ABC1MassMix", longname = "Absorbing Carbon Type 1 Mixing Ratio", stdname = None, vartype = "3d", ramsname = "ABC1MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "ABC2MassMix", longname = "Absorbing Carbon Type 2 Mixing Ratio", stdname = None, vartype = "3d", ramsname = "ABC2MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "SaltFilmMassMix", longname = "Sea Salt Film Drop Mixing Ratio", stdname = None, vartype = "3d", ramsname = "SALT_FILM_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "SaltJetMassMix", longname = "Sea Salt Jet Drop Mixing Ratio", stdname = None, vartype = "3d", ramsname = "SALT_JET_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "SaltSpumeMassMix", longname = "Sea Salt Spume Drop Mixing Ratio", stdname = None, vartype = "3d", ramsname = "SALT_SPUM_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "SubAeroRegenMassMix", longname = "Sub-Micron Regenerated Aerosol Mixing Ratio", stdname = None, vartype = "3d", ramsname = "REGEN_AERO1_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "SuperAeroRegenMassMix", longname = "Super-Micron Regenerated Aerosol Mixing Ratio", stdname = None, vartype = "3d", ramsname = "REGEN_AERO2_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "INPMassMix", longname = "Ice Nucleating Particle Mixing Ratio", stdname = None, vartype = "3d", ramsname = "RIFNP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))

#3D Aerosol Tracking Variables
    varlist.append(outvar(varname = "INPNucNumberMix", longname = "Number Mixing Ratio of INP Already Nucleated", stdname = None, vartype = "3d", ramsname = "IFNNUCP", unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "CloudINPNumberMix", longname = "Number Mixing Ratio of INP Within Cloud Droplets", stdname = None, vartype = "3d", ramsname = "IMMERCP", unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "RainINPNumberMix", longname = "Number Mixing Ratio of Within Rain", stdname = None, vartype = "3d", ramsname = "IMMERRP", unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    varlist.append(outvar(varname = "CloudAerosolMix", longname = "Mixing Ratio of Total Aerosol in Cloud Droplets", stdname = None, vartype = "3d", ramsname = "CNMCP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "DrizzleAerosolMix", longname = "Mixing Ratio of Total Aerosol in Drizzle", stdname = None, vartype = "3d", ramsname = "CNMDP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "RainAerosolMix", longname = "Mixing Ratio of Total Aerosol in Rain", stdname = None, vartype = "3d", ramsname = "CNMRP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "PrisAerosolMix", longname = "Mixing Ratio of Total Aerosol in Pristine Ice", stdname = None, vartype = "3d", ramsname = "CNMPP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "SnowAerosolMix", longname = "Mixing Ratio of Total Aerosol in Snow", stdname = None, vartype = "3d", ramsname = "CNMSP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "AggAerosolMix", longname = "Mixing Ratio of Total Aerosol in Aggregates", stdname = None, vartype = "3d", ramsname = "CNMAP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "GraupelAerosolMix", longname = "Mixing Ratio of Total Aerosol in Graupel", stdname = None, vartype = "3d", ramsname = "CNMGP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "HailAerosolMix", longname = "Mixing Ratio of Total Aerosol in Hail", stdname = None, vartype = "3d", ramsname = "CNMHP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "CloudDustMix", longname = "Mixing Ratio of Dust in Cloud Droplets", stdname = None, vartype = "3d", ramsname = "DNMCP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "DrizzleDustMix", longname = "Mixing Ratio of Dust in Drizzle", stdname = None, vartype = "3d", ramsname = "DNMDP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "RainDustMix", longname = "Mixing Ratio of Dust in Rain", stdname = None, vartype = "3d", ramsname = "DNMRP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "PrisDustMix", longname = "Mixing Ratio of Dust in Pristine Ice", stdname = None, vartype = "3d", ramsname = "DNMPP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "SnowDustMix", longname = "Mixing Ratio of Dust in Snow", stdname = None, vartype = "3d", ramsname = "DNMSP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "AggDustMix", longname = "Mixing Ratio of Dust in Aggregates", stdname = None, vartype = "3d", ramsname = "DNMAP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "GraupelDustMix", longname = "Mixing Ratio of Dust in Graupel", stdname = None, vartype = "3d", ramsname = "DNMGP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "HailDustMix", longname = "Mixing Ratio of Dust in Hail", stdname = None, vartype = "3d", ramsname = "DNMHP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "CloudIceNucDust", longname = "Mixing Ratio of Dust in Cloud Droplets from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINCP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "DrizzleIceNucDust", longname = "Mixing Ratio of Dust in Drizzle from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINDP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "RainIceNucDust", longname = "Mixing Ratio of Dust in Rain from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINRP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "PrisIceNucDust", longname = "Mixing Ratio of Dust in Pristine Ice from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINPP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "SnowIceNucDust", longname = "Mixing Ratio of Dust in Snow from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINSP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "AggIceNucDust", longname = "Mixing Ratio of Dust in Aggregates from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINAP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "GraupelIceNucDust", longname = "Mixing Ratio of Dust in Graupel from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINGP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "HailIceNucDust", longname = "Mixing Ratio of Dust in Hail from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINHP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "CloudSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Cloud Droplets", stdname = None, vartype = "3d", ramsname = "SNMCP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "DrizzleSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Drizzle", stdname = None, vartype = "3d", ramsname = "SNMDP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "RainSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Rain", stdname = None, vartype = "3d", ramsname = "SNMRP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "PrisSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Pristine Ice", stdname = None, vartype = "3d", ramsname = "SNMPP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "SnowSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Snow", stdname = None, vartype = "3d", ramsname = "SNMSP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "AggSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Aggregates", stdname = None, vartype = "3d", ramsname = "SNMAP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "GraupelSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Graupel", stdname = None, vartype = "3d", ramsname = "SNMGP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "HailSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Hail", stdname = None, vartype = "3d", ramsname = "SNMHP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "RegenSolubleSubAeroMix", longname = "Mixing Ratio of Sub-Micron Regenerated Soluble Aerosol", stdname = None, vartype = "3d", ramsname = "RESOL_AERO1_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    varlist.append(outvar(varname = "RegenSolubleSuperAeroMix", longname = "Mixing Ratio of Super-Micron Regenerated Soluble Aerosol", stdname = None, vartype = "3d", ramsname = "RESOL_AERO2_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))

#3D Precipitation Variables
    varlist.append(outvar(varname = "DrizzlePrecip3d", longname = "3D Rate of Drizzle Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVD", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "RainPrecip3d", longname = "3D Rate of Rain Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVR", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "PrisPrecip3d", longname = "3D Rate of Pristine Ice Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVP", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "SnowPrecip3d", longname = "3D Rate of Snow Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVS", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "AggPrecip3d", longname = "3D Rate of Aggregate Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVA", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "GraupelPrecip3d", longname = "3D Rate of Graupel Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVG", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "HailPrecip3d", longname = "3D Rate of Hail Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVH", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "PlatePrecip3d", longname = "3D Rate of Plate Ice Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVIP", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "ColumnIcePrecip3d", longname = "3D Rate of Column Ice Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVIC", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "DendritePrecip3d", longname = "3D Rate of Dendritic Snowflake Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVIP", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))

#Default Momentum Budget Variables
    varlist.append(outvar(varname = "ThetaBuoyDefault", longname = "Default Model Vertical Acceleration from Potential Temperature Difference from Reference Potential Temperature", stdname = None, vartype = "3d", ramsname = "WP_BUOY_THETA", invar = None, unitfactor = "defaultmomentumbudget", decnum = 3, units = "m s**-2"))
    varlist.append(outvar(varname = "CondBuoyDefault", longname = "Default Model Vertical Acceleration from Condensate Loading", stdname = None, vartype = "3d", ramsname = "WP_BUOY_COND", invar = None, unitfactor = "defaultmomentumbudget", decnum = 3, units = "m s**-2"))
    varlist.append(outvar(varname = "AdvDifDefault", longname = "Default Model Vertical Acceleration from Advection and Diffusion", stdname = None, vartype = "3d", ramsname = "WP_ADVDIF", invar = None, unitfactor = "defaultmomentumbudget", decnum = 3, units = "m s**-2"))

#Tracer Variables
    for i in range(1, 51):
        varlist.append(outvar(varname = f"Tracer{str(i).zfill(3)}Mix", longname = f"Number Mixing Ratio of Tracer Number {str(i).zfill(3)}", stdname = None, nptype = float64, vartype = "3d", ramsname = f"TRACERP{str(i).zfill(3)}", unitfactor = 10**-6, decnum = 6, units = "# mg**-1"))

#2D Grid Variables
    varlist.append(outvar(varname = "Lat2D", longname = "2D Latitude Field from Grid", stdname = "grid_latitude", vartype = "2d", ramsname = "GLAT", invar = None, unitfactor = 1, decnum = 3, units = "degree"))
    varlist.append(outvar(varname = "Lon2D", longname = "2D Longitude Field from Grid", stdname = "grid_longitude", vartype = "2d", ramsname = "GLON", invar = None, unitfactor = 1, decnum = 3, units = "degree"))
    varlist.append(outvar(varname = "Topo", longname = "Topographic Height Above Mean Sea Level", stdname = "surface_altitude", vartype = "2d", ramsname = "TOPT", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    varlist.append(outvar(varname = "TopoRoughness", longname = "Topographic Roughness Length", stdname = None, vartype = "2d", ramsname = "TOPZO", invar = None, unitfactor = 1, decnum = 3, units = "m"))

#2D Surface File Input Variables
    varlist.append(outvar(varname = "SoilType", longname = "Soil Type of Current Soil Level", stdname = "soil_type", vartype = "2dSoil", ramsname = "SOIL_TEXT", invar = None, unitfactor = 1, decnum = 0, units = "Unitless"))
    varlist.append(outvar(varname = "VegClass", longname = "Vegetation Type over Current Patch", stdname = "land_cover", vartype = "2dPatch", ramsname = "LEAF_CLASS", invar = None, unitfactor = 1, decnum = 0, units = "Unitless"))
    varlist.append(outvar(varname = "PatchArea", longname = "Fraction of Grid Cell Occupied by Current Patch", stdname = None, vartype = "2dPatch", ramsname = "PATCH_AREA", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "PastNDVI", longname = "Past Vegetation Normalized Difference Vegetation Index over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "VEG_NDVIP", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "NDVI", longname = "Vegetation Normalized Difference Vegetation Index over Current Patch", stdname = "normalized_difference_vegetation_index", vartype = "2dPatch", ramsname = "VEG_NDVIC", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "FutureNDVI", longname = "Future Vegetation Normalized Difference Vegetation Index over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "VEG_NDVIF", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))

#2D Surface Characteristics 
    varlist.append(outvar(varname = "PatchRough", longname = "Patch Roughness Length", stdname = None, vartype = "2dPatch", ramsname = "PATCH_ROUGH", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    varlist.append(outvar(varname = "SoilRough", longname = "Soil Type Contribution to Roughness Length over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "SOIL_ROUGH", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    varlist.append(outvar(varname = "VegRough", longname = "Vegetation Type Contribution to Roughness Length over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "VEG_ROUGH", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    varlist.append(outvar(varname = "VegFrac", longname = "Fraction of Current Patch Covered by Vegetation", stdname = "vegetation_area_fraction", vartype = "2dPatch", ramsname = "VEG_FRACAREA", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    varlist.append(outvar(varname = "VegLeafArea", longname = "Green Leaf Area Index (m^2 of Green Leaves to m^2 of Surface)", stdname = "leaf_area_index", vartype = "2dPatch", ramsname = "VEG_LAI", invar = None, unitfactor = 1, decnum = 3, units = "m**2 m**-2"))
    varlist.append(outvar(varname = "VegTotalArea", longname = "Total Vegetation Area Index (m^2 of Green Leaves+Dead Leaves+Stems+ Trunks to m^2 of Surface)", stdname = None, vartype = "2dPatch", ramsname = "VEG_TAI", invar = None, unitfactor = 1, decnum = 3, units = "m**2 m**-2"))
    varlist.append(outvar(varname = "VegHeight", longname = "Vegetation Height Above Ground", stdname = None, vartype = "2dPatch", ramsname = "VEG_HEIGHT", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    varlist.append(outvar(varname = "VegAlbedo", longname = "Vegetation Broadband Albedo", stdname = None, vartype = "2dPatch", ramsname = "VEG_ALBEDO", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "SoilMoisture", longname = "Volumetric Soil Moisture in Current Soil Level", stdname = "volume_fraction_of_condensed_water_in_soil", vartype = "2dSoil", ramsname = "SOIL_WATER", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "SoilEnergy", longname = "Soil Energy With Respect to Freezing Within Current Soil Level", stdname = None, vartype = "2dSoil", ramsname = "SOIL_ENERGY", invar = None, unitfactor = 1, decnum = 3, units = "J kg**-1"))
    varlist.append(outvar(varname = "SnowLevelCount", longname = "Number of Snow Levels over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "SFCWATER_NLEV", invar = None, unitfactor = 1, decnum = 0, units = "Unitless"))
    varlist.append(outvar(varname = "SnowLevelMass", longname = "Snow Mass Within Current Snow Level", stdname = None, vartype = "2dSnow", ramsname = "SFCWATER_MASS", invar = None, unitfactor = 1, decnum = 3, units = "kg m**-2"))
    varlist.append(outvar(varname = "SnowLevelEnergy", longname = "Snow Energy With Respect to Freezing Within Current Snow Level", stdname = None, vartype = "2dSnow", invar = None, ramsname = "SFCWATER_ENERGY", unitfactor = 1, decnum = 3, units = "J m**-3"))
    varlist.append(outvar(varname = "SnowLevelDepth", longname = "Snow Depth Within Current Snow Level", stdname = None, vartype = "2dSnow", ramsname = "SFCWATER_DEPTH", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    varlist.append(outvar(varname = "UStar", longname = "Friction Velocity for Surface Fluxes over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "USTAR", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))
    varlist.append(outvar(varname = "TStar", longname = "Scaling Temperature for Surface Fluxes over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "TSTAR", invar = None, unitfactor = 1, decnum = 3, units = "K"))
    varlist.append(outvar(varname = "RStar", longname = "Scaling Vapor Mixing Ratio for Surface Fluxes over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "RSTAR", invar = None, unitfactor = 1000, decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "StomatalResistance", longname = "Vegetation Stomatal Resistance", stdname = None, vartype = "2dPatch", ramsname = "STOM_RESIST", invar = None, unitfactor = 1, decnum = 2, units = "s m**-1"))
    varlist.append(outvar(varname = "VegSurfaceWater", longname = "Liquid Water Depth on Vegetation Surface", stdname = "canopy_water_amount", vartype = "2dPatch", ramsname = "VEG_WATER", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    varlist.append(outvar(varname = "VegTemp", longname = "Vegetation Temperature", stdname = "canopy_temperature", vartype = "2dPatch", ramsname = "VEG_TEMP", invar = None, unitfactor = 1, decnum = 3, units = "K"))
    varlist.append(outvar(varname = "CanopyVaporMix", longname = "Vapor Mixing Ratio in Vegetation Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "CAN_RVAP", invar = None, unitfactor = 1000, decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "CanopyTemp", longname = "Temperature in Vegetation Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "CAN_TEMP", invar = None, unitfactor = 1, decnum = 3, units = "K"))
    varlist.append(outvar(varname = "GroundSatVaporMix", longname = "Saturation Vapor Mixing Ratio of Soil/Snow Surface (Used for Snow Sublimation and Dew Formation)", stdname = None, vartype = "2dPatch", ramsname = "GROUND_RSAT", invar = None, unitfactor = 1000, decnum = 3, units = "g kg**-1"))
    varlist.append(outvar(varname = "SoilVaporMix", longname = "Vapor Mixing Ratio of Soil Surface (Used for Evaporation from Soil When No Snow Present)", stdname = None, vartype = "2dPatch", ramsname = "GROUND_RVAP", invar = None, unitfactor = 1000, decnum = 3, units = "g kg**-1"))

#2D SiB Variables
    varlist.append(outvar(varname = "CO2Mix", longname = "CO2 Mass Fraction (Divide by 1.5172*10**-6 to get CO2 ppm)", stdname = "mass_fraction_of_carbon_dioxide_in_air", vartype = "2dPatch", ramsname = "RCO2P", invar = None, unitfactor = 1, decnum = 7, units = "kg kg**-1"))
    varlist.append(outvar(varname = "VegSnow", longname = "Mass of Snow Covering Vegetation", stdname = "canopy_snow_amount", vartype = "2dPatch", ramsname = "SNOW1", invar = None, unitfactor = 1, decnum = 4, units = "kg m**-2"))
    varlist.append(outvar(varname = "SoilSnow", longname = "Mass of Snow Covering Bare Ground", stdname = "surface_snow_amount", vartype = "2dPatch", ramsname = "SNOW2", invar = None, unitfactor = 1, decnum = 4, units = "kg m**-2"))
    varlist.append(outvar(varname = "VegStore", longname = "Vegetation Liquid Store", stdname = None, vartype = "2dPatch", ramsname = "CAPAC1", invar = None, unitfactor = 1, decnum = 4, units = "kg m**-2"))
    varlist.append(outvar(varname = "SoilStore", longname = "Soil Surface Liquid Store", stdname = None, vartype = "2dPatch", ramsname = "CAPAC2", invar = None, unitfactor = 1, decnum = 4, units = "kg m**-2"))
    varlist.append(outvar(varname = "CanopyCO2", longname = "Canopy Air Space CO2 Partial Pressure", stdname = None, vartype = "2dPatch", ramsname = "PCO2AP", invar = None, unitfactor = 1, decnum = 2, units = "Pa"))
    varlist.append(outvar(varname = "CO2Flux", longname = "CO2 Flux from Canopy Air Space to Atmosphere", stdname = None, vartype = "2dPatch", ramsname = "CO2FLX", invar = None, unitfactor = 1, decnum = 4, units = "mol m**-2 s**-1"))
    varlist.append(outvar(varname = "SIBSurfaceAlbedo", longname = "Broadband Surface Albedo of Current Patch", stdname = None, vartype = "2dPatch", ramsname = "SFCSWA", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "SIBLongwaveUpsrf", longname = "Upwelling Longwave Radiation From the Surface over the Current Patch", stdname = None, vartype = "2dPatch", ramsname = "UPLWRF", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "CO2Uptake", longname = "Uptake of CO2 by Canopy Plants", stdname = None, vartype = "2dPatch", ramsname = "ASSIMN", invar = None, unitfactor = 1, decnum = 3, units = "umol m**-2 s**-1"))
    varlist.append(outvar(varname = "GroundRespiration", longname = "Ground Respiration Flux", stdname = None, vartype = "2dPatch", ramsname = "RESPGF", invar = None, unitfactor = 1, decnum = 3, units = "umol m**-2 s**-1"))
    varlist.append(outvar(varname = "LeafResistStress", longname = "Leaf Surface Humidity Resistance Stress", stdname = None, vartype = "2dPatch", ramsname = "RSTFAC1", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "SoilResistStress", longname = "Soil Moisture Resistance Stress", stdname = None, vartype = "2dPatch", ramsname = "RSTFAC2", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "TempResistStress", longname = "Temperature Resistance Stress", stdname = None, vartype = "2dPatch", ramsname = "RSTFAC3", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "Transpiration", longname = "Latent Heat Flux into Canopy Air Space from Transpiration", stdname = "transpiration_flux", vartype = "2dPatch", ramsname = "ECT", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "CanopyIntercept", longname = "Latent Heat Flux into Canopy Air Space from Evaporation of Water on Vegetation Surface", stdname = "tendency_of_canopy_water_amount_from_evaporation_of_intercepted_precipitation", vartype = "2dPatch", ramsname = "ECI", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "GroundIntercept", longname = "Latent Heat Flux into Canopy Air Space from Evaporation of Water on Soil Surface", stdname = None, vartype = "2dPatch", ramsname = "EGI", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "GroundEvap", longname = "Latent Heat Flux into Canopy Air Space from Evaporation of Soil Moisture", stdname = "water_evaporation_flux_from_soil", vartype = "2dPatch", ramsname = "EGS", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "VegSensHeatFlux", longname = "Sensible Heat Flux into Canopy Air Space from Vegetation", stdname = None, vartype = "2dPatch", ramsname = "HC", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "SoilSensHeatFlux", longname = "Sensible Heat Flux into Canopy Air Space from Soil Surface", stdname = None, vartype = "2dPatch", ramsname = "HG", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "CanopyToAtmosResistance", longname = "Aerodynamic Resistance from Canopy Air Space to Atmosphere", stdname = None, vartype = "2dPatch", ramsname = "RA", invar = None, unitfactor = 1, decnum = 2, units = "s m**-1"))
    varlist.append(outvar(varname = "VegToCanopyResistance", longname = "Aerodynamic Resistance from Vegetation Surface to Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "RB", invar = None, unitfactor = 1, decnum = 2, units = "s m**-1"))
    varlist.append(outvar(varname = "TotalCanopyResistance", longname = "Total Aerodynamic Resistance of Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "RC", invar = None, unitfactor = 1, decnum = 2, units = "s m**-1"))
    varlist.append(outvar(varname = "SoilToCanopyResistance", longname = "Aerodynamic Resistance from Soil Surface to Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "RD", invar = None, unitfactor = 1, decnum = 2, units = "s m**-1"))
    varlist.append(outvar(varname = "Runoff", longname = "Depth of Surface and Subsurface Runoff", stdname = None, vartype = "2dPatch", ramsname = "ROFF", invar = None, unitfactor = 1, decnum = 2, units = "mm"))
    varlist.append(outvar(varname = "GreenFrac", longname = "Greenness Fraction of Vegetation (Area of Green Leaves to Total Plant Area)", stdname = None, vartype = "2dPatch", ramsname = "GREEN", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "PhotoAbsorb", longname = "Fraction of Photosynthetically Active Radiation Absorbed by Vegetation", stdname = None, vartype = "2dPatch", ramsname = "APAR", invar = None, unitfactor = 1, decnum = 4, units = "Unitless"))
    varlist.append(outvar(varname = "VentilationMassFlux", longname = "Ventilation Mass Flux", stdname = None, vartype = "2dPatch", ramsname = "VENTMF", invar = None, unitfactor = 1, decnum = 4, units = "kg m**-2 s**-1"))
    varlist.append(outvar(varname = "ChloroCO2Pressure", longname = "Partial Pressure of CO2 in Chloroplasts", stdname = None, vartype = "2dPatch", ramsname = "PCO2C", invar = None, unitfactor = 1, decnum = 2, units = "Pa"))
    varlist.append(outvar(varname = "LeafInternalCO2Pressure", longname = "Partial Pressure of CO2 in Internal Leaf Space", stdname = None, vartype = "2dPatch", ramsname = "PCO2I", invar = None, unitfactor = 1, decnum = 2, units = "Pa"))
    varlist.append(outvar(varname = "LeafSurfaceCO2Pressure", longname = "Partial Pressure of CO2 at Leaf Surface", stdname = None, vartype = "2dPatch", ramsname = "PCO2S", invar = None, unitfactor = 1, decnum = 2, units = "Pa"))
    varlist.append(outvar(varname = "CanopyVaporPres", longname = "Water Vapor Pressure in Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "EA", invar = None, unitfactor = 1, decnum = 3, units = "hPa"))
    varlist.append(outvar(varname = "ReferenceVaporPres", longname = "Reference Level Water Vapor Pressure", stdname = None, vartype = "2dPatch", ramsname = "EM", invar = None, unitfactor = 1, decnum = 3, units = "hPa"))
    varlist.append(outvar(varname = "CanopyHumidity", longname = "Relative Humidity of Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "RHA", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "VisibleDiffuse", longname = "Visible Diffuse Radiation", stdname = "surface_diffuse_downwelling_shortwave_flux_in_air", vartype = "2dPatch", ramsname = "RADVBC", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "VisibleDirect", longname = "Visible Direct Radiation", stdname = "surface_direct_downwelling_shortwave_flux_in_air", vartype = "2dPatch", ramsname = "RADVDC", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "NIRDiffuse", longname = "Near-Infrared Diffuse Radiation", stdname = None, vartype = "2dPatch", ramsname = "RADNBC", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "NIRDirect", longname = "Near-Infrared Direct Radiation", stdname = None, vartype = "2dPatch", ramsname = "RADNDC", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "PsyConstant", longname = "Psychrometric Constant", stdname = None, vartype = "2dPatch", ramsname = "PSY", invar = None, unitfactor = 1, decnum = 3, units = "hPa deg**-1"))

#2D Precipitation Variables
    varlist.append(outvar(varname = "DrizzlePrecipRate", longname = "Rate of Drizzle Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRD", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "RainPrecipRate", longname = "Rate of Rain Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = "rainfall_rate", vartype = "2d", ramsname = "PCPRR", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "PrisPrecipRate", longname = "Rate of Pristine Ice Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRP", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "SnowPrecipRate", longname = "Rate of Snow Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRS", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "AggPrecipRate", longname = "Rate of Aggregate Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRA", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "GraupelPrecipRate", longname = "Rate of Graupel Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRG", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "HailPrecipRate", longname = "Rate of Hail Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRH", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "PlatePrecipRate", longname = "Rate of Plate Ice Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRIP", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "ColumnIcePrecipRate", longname = "Rate of Column Ice Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRIC", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "DendritePrecipRate", longname = "Rate of Dendritic Snowflake Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRID", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "DrizzlePrecipTotal", longname = "Total Accumulated Drizzle Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPD", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    varlist.append(outvar(varname = "RainPrecipTotal", longname = "Total Accumulated Rain Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPR", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    varlist.append(outvar(varname = "PrisPrecipTotal", longname = "Total Accumulated Pristine Ice Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPP", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    varlist.append(outvar(varname = "SnowPrecipTotal", longname = "Total Accumulated Snow Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPS", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    varlist.append(outvar(varname = "AggPrecipTotal", longname = "Total Accumulated Aggregate Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPA", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    varlist.append(outvar(varname = "GraupelPrecipTotal", longname = "Total Accumulated Graupel Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPG", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    varlist.append(outvar(varname = "HailPrecipTotal", longname = "Total Accumulated Hail Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPH", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    varlist.append(outvar(varname = "PlatePrecipTotal", longname = "Total Accumulated Plate Ice Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPIP", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    varlist.append(outvar(varname = "ColumnIcePrecipTotal", longname = "Total Accumulated Column Ice Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPIC", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    varlist.append(outvar(varname = "DendritePrecipTotal", longname = "Total Accumulated Dendritic Snowflake Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPID", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    varlist.append(outvar(varname = "LEAFPrecipRate", longname = "Instantaneous Mass Flux of Microphysics Precipitation (used by LEAF/SIB)", stdname = None, vartype = "2d", ramsname = "PCPG", invar = None, unitfactor = "pertimestep", decnum = 3, units = "kg m**-2"))
    varlist.append(outvar(varname = "LEAFPrecipEnergy", longname = "Instantaneous Energy Flux of Microphysics Precipitation (used by LEAF/SIB)", stdname = None, vartype = "2d", ramsname = "QPCPG", invar = None, unitfactor = "pertimestep", decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "LEAFPrecipDepth", longname = "Instantaneous Rate of Precipitation Depth (used by LEAF/SIB)", stdname = None, vartype = "2d", ramsname = "DPCPG", invar = None, unitfactor = "pertimestep", decnum = 3, units = "m s**-1"))

#Aerosol Deposition
    varlist.append(outvar(varname = "TotalDustDeposition", longname = "Total Deposition of Dust on the Surface", stdname = None, vartype = "2d", ramsname = "ACCPDUST", invar = None, unitfactor = 1, decnum = 3, units = "kg m**-2"))
    varlist.append(outvar(varname = "TotalAerosolDeposition", longname = "Total Deposition of Aerosol on the Surface", stdname = None, vartype = "2d", ramsname = "ACCPAERO", invar = None, unitfactor = 1, decnum = 3, units = "kg m**-2"))
    varlist.append(outvar(varname = "DustDepositionRate", longname = "Rate of Deposition of Dust on the Surface", stdname = None, vartype = "2d", ramsname = "PCPRDUST", invar = None, unitfactor = 1, decnum = 3, units = "kg m**-2 s**-1"))
    varlist.append(outvar(varname = "AeroDepositionRate", longname = "Rate of Deposition of Aerosol on the Surface", stdname = None, vartype = "2d", ramsname = "PCPRAERO", invar = None, unitfactor = 1, decnum = 3, units = "kg m**-2"))

#2D Radiation Variables
    varlist.append(outvar(varname = "ShortwaveDownSrf", longname = "Surface Downwelling Shortwave Radiation", stdname = "surface_downwelling_shortwave_flux", vartype = "2d", ramsname = "RSHORT", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "LongwaveDownSrf", longname = "Surface Downwelling Longwave Radiation", stdname = "surface_downwelling_longwave_flux", vartype = "2d", ramsname = "RLONG", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "LongwaveUpSrf", longname = "Surface Upwelling Longwave Radiation", stdname = "surface_upwelling_longwave_flux", vartype = "2d", ramsname = "RLONGUP", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "AerosolOpticalDepth", longname = "Aerosol Optical Depth in Visible Light", stdname = "atmosphere_optical_thickness_due_to_ambient_aerosol_particles", vartype = "2d", ramsname = "AODT", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "SurfaceAlbedo", longname = "Surface Broadband Albedo", stdname = "surface_albedo", vartype = "2d", ramsname = "ALBEDT", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    varlist.append(outvar(varname = "SolarZenithAngle", longname = "Solar Zenith Angle", stdname = "surface_albedo", vartype = "2d", ramsname = "COSZ", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))

#2D Surface Fluxes
    varlist.append(outvar(varname = "SensibleHeatFlux", longname = "Surface Sensible Heat Flux", stdname = "surface_upward_sensible_heat_flux", vartype = "2d", ramsname = "SFLUX_T", invar = None, unitfactor = 1004, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "LatentHeatFlux", longname = "Surface Latent Heat Flux", stdname = "surface_upward_latent_heat_flux", vartype = "2d", ramsname = "SFLUX_R", invar = None, unitfactor = 2.5*10**6, decnum = 3, units = "W m**-2"))
    varlist.append(outvar(varname = "UMomentumFlux", longname = "Surface U-Momentum Flux", stdname = None, vartype = "2d", ramsname = "SFLUX_U", invar = None, unitfactor = 1, decnum = 3, units = "Pa"))
    varlist.append(outvar(varname = "VMomentumFlux", longname = "Surface V-Momentum Flux", stdname = None, vartype = "2d", ramsname = "SFLUX_V", invar = None, unitfactor = 1, decnum = 3, units = "Pa"))
    varlist.append(outvar(varname = "WMomentumFlux", longname = "Surface W-Momentum Flux", stdname = None, vartype = "2d", ramsname = "SFLUX_W", invar = None, unitfactor = 1, decnum = 3, units = "Pa"))

#2D Convective Parameterization Variables
    varlist.append(outvar(varname = "ConvectivePrecipTotal", longname = "Total Accumulated Precipitation at the Surface from Convective Parameterization", stdname = "convective_precipitation_amount", vartype = "2d", ramsname = "ACONPR", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    varlist.append(outvar(varname = "ConvectivePrecipRate", longname = "Rate of Precipitation at the Surface from Convective Parameterization", stdname = "convective_precipitation_rate", vartype = "2d", ramsname = "CONPRR", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    varlist.append(outvar(varname = "ConvectiveTimesteps", longname = "Number of Timesteps at Horizontal Gridpoint Maintaining Convective Behavior", stdname = None, vartype = "2d", ramsname = "NCA", invar = None, unitfactor = 1, decnum = 3, units = 1))
    varlist.append(outvar(varname = "ConvectiveInitCheck", longname = "Check if Pre-Convection Checks at Horizontal Gridpoint Satisfied", stdname = None, vartype = "2d", ramsname = "CONVGO", invar = None, unitfactor = 1, decnum = 3, units = 1))
    
    
    return varlist

# def datavarinit(varentry, ftime, coords):
#     pass

def gen_vardict(uservarlist, datavarlist, fullvarlist):
    vdict = {}
    uservarset = set(uservarlist)
    datavarset = set(datavarlist)
    for entry in fullvarlist:
        if entry.ramsname in (uservarset & datavarset):
            vdict[entry.ramsname] = entry
        elif entry.ramsname in uservarset and entry.ramsname not in datavarset:
            print(f"User requested variable {entry.ramsname}, but {entry.ramsname} is not present in the RAMS output file!")
    return vdict