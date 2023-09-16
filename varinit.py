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
        #stdname (string) is the "standard name" according to the cf naming conventions, to be used in the xarray attributes
        #vector is whether this variable is a vector (u,v,w) that needs a special interpolation routine.
        #offline (bool) is whether this variable is simply a renamed and truncated version of a raw RAMS variable (ex. Theta, SrfTemp) or an offline-calculated variable (like pressure or cloud top height)
        #vartype (string) is whether this variable is 3d (z,y,x), 2d (y,x), 2d-patch (patch, y, x), 2d-snow (patch, snowlev, y, x), 2d-soil (patch, soillev, y, x)
        #ntype (numpy dtype) is the dtype to use (float32 or float64) for the variable
        #ramsname (string) is the name of the raw rams variable read in before truncation and renaming. Only applies to online variables (Theta, VaporMix)
        #invar (list of strings) is the other variables (if any) that need to be passed in to calculate this variable. Only applies for offline-calculated variables.
        #decnum (int) is the number of decimal points to use for storing this variable. This will only come into play after all unit conversions
        #unitfactor (float or string) is the factor by which this variable needs to be modified by its original units (for example, 2.5*10^6 for SFLUX_R (kg/m^2*s) to LatHeatFlux (W/m^2) or 1000 for RCP (kg/kg) to CloudMassMix (g/kg), 1000/afiletime for AGGREGATET to Aggregation (change in mixing ratio due to aggregation in the previous timestep))
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

def fillvarlist(uservars, datavarlist):
    #fillvarlist is a somewhat expensive function, dealing with lots of classes. So, this is just run at the start of the post-processing (NOT for every timestep). The post-processing routine assumes that you want the same variables for all times in the post-processing. If you want to analyze new variables at certain times, restart the program with a new starttime and a new variable list file.
    varlist = []

    #3D Vector variables
    if "UC" in uservars and "UC" in datavarlist:
        varlist.append(outvar(varname = "u", longname = "Eastward Wind Velocity", stdname = "eastward_wind", vector = True, vartype = "3d", nptype = float32, ramsname = "UC", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))
    if "VC" in uservars and "VC" in datavarlist:
        varlist.append(outvar(varname = "v", longname = "Northward Wind Velocity", stdname = "northward_wind", vector = True, vartype = "3d", nptype = float32, ramsname = "VC", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))
    if "WC" in uservars and "WC" in datavarlist:
        varlist.append(outvar(varname = "w", longname = "Vertical Velocity", stdname = "upward_air_velocity", vector = True, vartype = "3d", nptype = float32, ramsname = "WC", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))
    if "UP" in uservars and "UP" in datavarlist:
        varlist.append(outvar(varname = "PastU", longname = "Past Timestep Eastward Wind Velocity", stdname = None, vector = True, vartype = "3d", nptype = float32, ramsname = "UP", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))
    if "VP" in uservars and "VP" in datavarlist:
        varlist.append(outvar(varname = "PastV", longname = "Past Timestep Northward Wind Velocity", stdname = None, vector = True, vartype = "3d", nptype = float32, ramsname = "VP", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))
    if "WP" in uservars and "WP" in datavarlist:
        varlist.append(outvar(varname = "PastW", longname = "Past Timestep Vertical Velocity", stdname = None, vector = True, vartype = "3d", nptype = float32, ramsname = "WP", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))

    #3D Meteorology Scalar Variables
    if "RV" in uservars and "RV" in datavarlist:
        varlist.append(outvar(varname = "VaporMix", longname = "Vapor Mixing Ratio", stdname = "humidity_mixing_ratio", vartype = "3d", nptype = float32, ramsname = "RV", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    if "RTP" in uservars and "RTP" in datavarlist:
        varlist.append(outvar(varname = "WaterMix", longname = "Water Mixing Ratio", stdname = "water_mixing_ratio", vartype = "3d", nptype = float32, ramsname = "RTP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    if "THETA" in uservars and "THETA" in datavarlist:
        varlist.append(outvar(varname = "Theta", longname = "Potential Temperature", stdname = "air_potential_temperature", vartype = "3d", nptype = float32, ramsname = "THETA", invar = None, unitfactor = 1, decnum = 3, units = "K"))
    if "PI" in uservars and "PI" in datavarlist:
        varlist.append(outvar(varname = "Exner", longname = "Exner Function * Cp", stdname = None, vartype = "3d", nptype = float32, ramsname = "PI", invar = None, unitfactor = 1, decnum = 3, units = 1))
    if "DN0" in uservars and "DN0" in datavarlist:
        varlist.append(outvar(varname = "ReferenceDensity", longname = "Reference Air Density", stdname = None, vartype = "3d", nptype = float32, ramsname = "DN0", invar = None, unitfactor = 1, decnum = 3, units = "kg m**-3"))
    if "PI0" in uservars and "PI0" in datavarlist:
        varlist.append(outvar(varname = "ReferenceExner", longname = "Reference Exner Function * Cp", stdname = None, vartype = "3d", nptype = float32, ramsname = "PI0", invar = None, unitfactor = 1, decnum = 3, units = 1))
    if "TH0" in uservars and "TH0" in datavarlist:
        varlist.append(outvar(varname = "ReferenceTheta", longname = "Reference Potential Temperature", stdname = None, vartype = "3d", nptype = float32, ramsname = "TH0", invar = None, unitfactor = 1, decnum = 3, units = "K"))
    if "PP" in uservars and "PP" in datavarlist:
        varlist.append(outvar(varname = "PastExnerPrime", longname = "Past Perturbation Exner Function * Cp", stdname = None, vartype = "3d", ramsname = "PP", invar = None, unitfactor = 1, decnum = 3, units = 1))
    if "PC" in uservars and "PC" in datavarlist:
        varlist.append(outvar(varname = "ExnerPrime", longname = "Current Perturbation Exner Function * Cp", stdname = None, vartype = "3d", ramsname = "PC", invar = None, unitfactor = 1, decnum = 3, units = 1))
    if "THP" in uservars and "THP" in datavarlist:
        varlist.append(outvar(varname = "Thetail", longname = "Ice-Liquid Potential Temperature", stdname = None, vartype = "3d", ramsname = "THP", invar = None, unitfactor = 1, decnum = 3, units = "K"))

    #3D Turbulence
    if "TKEP" in uservars and "TKEP" in datavarlist:
        varlist.append(outvar(varname = "TKE", longname = "Turbulent Kinetic Energy", stdname = "specific_turbulent_kinetic_energy_of_air", vartype = "3d", ramsname = "TKEP", invar = None, unitfactor = 1, decnum = 3, units = "J kg**-1"))
    if "HKH" in uservars and "HKH" in datavarlist:
        varlist.append(outvar(varname = "HorizHeatDiffusivityUnWeighted", longname = "Horizontal Eddy Diffusivity Coefficient for Heat (No Density Weighting)", stdname = None, vartype = "3d", ramsname = "HKH", invar = None, unitfactor = 1, decnum = 3, units = "m**2 s**-1"))
    if "VKH" in uservars and "VKH" in datavarlist:
        varlist.append(outvar(varname = "VertHeatDiffusivityUnWeighted", longname = "Vertical Eddy Diffusivity Coefficient for Heat (No Density Weighting)", stdname = None, vartype = "3d", ramsname = "VKH", invar = None, unitfactor = 1, decnum = 3, units = "m**2 s**-1"))
    if "RVKH" in uservars and "RVKH" in datavarlist:
        varlist.append(outvar(varname = "VertHeatDiffusivityWeighted", longname = "Vertical Eddy Diffusivity Coefficient for Heat", stdname = "atmosphere_heat_diffusivity", vartype = "3d", ramsname = "RVKH", invar = None, unitfactor = 1, decnum = 3, units = "m**2 s**-1"))
    if "RHKM" in uservars and "RHKM" in datavarlist:
        varlist.append(outvar(varname = "HorizMomentumDiffusivityWeighted", longname = "Horizontal Eddy Diffusivity Coefficient for Momentum", stdname = None, vartype = "3d", ramsname = "RHKM", invar = None, unitfactor = 1, decnum = 3, units = "m**2 s**-1"))
    if "RVKM" in uservars and "RVKM" in datavarlist:
        varlist.append(outvar(varname = "VertMomentumDiffusivityWeighted", longname = "Vertical Eddy Diffusivity Coefficient for Momentum", stdname = "atmosphere_momentum_diffusivity", vartype = "3d", ramsname = "RVKM", invar = None, unitfactor = 1, decnum = 3, units = "m**2 s**-1"))
    
    #3d Convective Parameterization Variables
    if "THSRC" in uservars and "THSRC" in datavarlist:
        varlist.append(outvar(varname = "ConvectiveHeating", longname = "Rate of Change in Air Temperature from Convective Parameterization", stdname = "tendency_of_air_temperature_due_to_convection", vartype = "3d", ramsname = "THSRC", invar = None, unitfactor = 1, decnum = 4, units = "K s**-1"))
    if "RTSRC" in uservars and "RTSRC" in datavarlist:
        varlist.append(outvar(varname = "ConvectiveMoistening", longname = "Rate of Change in Vapor Mixing Ratio from Convective Parameterization", stdname = None, vartype = "3d", ramsname = "RTSRC", unitfactor = "gperkg", decnum = 4, units = "g kg**-1 s**-1"))
    if "RCSRC" in uservars and "RCSRC" in datavarlist:
        varlist.append(outvar(varname = "ConvectiveCloudMixChange", longname = "Rate of Change in Cloud Droplet Mixing Ratio from Convective Parameterization", stdname = None, vartype = "3d", ramsname = "RCSRC", unitfactor = 10**6, decnum = 4, units = "mg kg**-1 s**-1"))
    if "RRSRC" in uservars and "RRSRC" in datavarlist:
        varlist.append(outvar(varname = "ConvectiveRainMixChange", longname = "Rate of Change in Rain Mixing Ratio from Convective Parameterization", stdname = None, vartype = "3d", ramsname = "RRSRC", unitfactor = 10**6, decnum = 4, units = "mg kg**-1 s**-1"))
    if "RPSRC" in uservars and "RPSRC" in datavarlist:
        varlist.append(outvar(varname = "ConvectivePrisMixChange", longname = "Rate of Change in Pristine Ice Mixing Ratio from Convective Parameterization", stdname = None, vartype = "3d", ramsname = "RPSRC", unitfactor = 10**6, decnum = 4, units = "mg kg**-1 s**-1"))
    if "RSSRC" in uservars and "RSSRC" in datavarlist:
        varlist.append(outvar(varname = "ConvectiveSnowMixChange", longname = "Rate of Change in Snow Mixing Ratio from Convective Parameterization", stdname = None, vartype = "3d", ramsname = "RSSRC", unitfactor = 10**6, decnum = 4, units = "mg kg**-1 s**-1"))
    if "W0AVG" in uservars and "W0AVG" in datavarlist:
        varlist.append(outvar(varname = "ConvectiveVertVelocity", longname = "Running Mean Average of Vertical Velocity", stdname = None, vartype = "3d", ramsname = "W0AVG", unitfactor = 1, decnum = 3, units = "m s**-1"))
    if "W0AVGLT" in uservars and "W0AVGLT" in datavarlist:
        varlist.append(outvar(varname = "ConvectiveVertVelocityContra", longname = "Running Mean Average of Horizontal Components of Contravariant Vertical Velocity", stdname = None, vartype = "3d", ramsname = "W0AVGLT", unitfactor = 1, decnum = 3, units = "m s**-1"))

    #3D Radiation Variables
    if "BEXT" in uservars and "BEXT" in datavarlist:
        varlist.append(outvar(varname = "VisLength", longname = "Length of Visibility", stdname = "visibility_in_air", vartype = "3d", nptype = float32, ramsname = "BEXT", invar = None, unitfactor = 1000, decnum = 3, units = "m"))
    if "FTHRD" in uservars and "FTHRD" in datavarlist:
        varlist.append(outvar(varname = "RadiativeHeating", longname = "Radiative Heating Rate in Air", stdname = "tendency_of_air_temperature_due_to_radiative_heating", vartype = "3d", nptype = float32, ramsname = "FTHRD", invar = None, unitfactor = 1, decnum = 3, units = "K s**-1"))
    if "LWDN" in uservars and "LWDN" in datavarlist:
        varlist.append(outvar(varname = "LongwaveDown", longname = "Downwelling Longwave Radiative Flux", stdname = "downwelling_longwave_flux_in_air", vartype = "3d", nptype = float32, ramsname = "LWDN", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "LWUP" in uservars and "LWUP" in datavarlist:
        varlist.append(outvar(varname = "LongwaveUp", longname = "Upwelling Longwave Radiative Flux", stdname = "upwelling_longwave_flux_in_air", vartype = "3d", nptype = float32, ramsname = "LWUP", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "SWDN" in uservars and "SWUP" in datavarlist:
        varlist.append(outvar(varname = "ShortwaveDown", longname = "Downwelling Shortwave Radiative Flux", stdname = "downwelling_shortwave_flux_in_air", vartype = "3d", nptype = float32, ramsname = "SWDN", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "SWUP" in uservars and "SWUP" in datavarlist:
        varlist.append(outvar(varname = "ShortwaveUp", longname = "Upwelling Shortave Radiative Flux", stdname = "upwelling_shortwave_flux_in_air", vartype = "3d", nptype = float32, ramsname = "SWUP", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))

    #Hydrometeor Mass/Number Mixing Ratios
    if "CCP" in uservars and "CCP" in datavarlist:
        varlist.append(outvar(varname = "CloudNumberMix", longname = "Cloud Droplet Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CCP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    if "CDP" in uservars and "CDP" in datavarlist:
        varlist.append(outvar(varname = "DrizzleNumberMix", longname = "Drizzle Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CDP", invar = None, unitfactor = 1, decnum = 0, units = "# kg**-1"))
    if "CRP" in uservars and "CRP" in datavarlist:
        varlist.append(outvar(varname = "RainNumberMix", longname = "Rain Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CRP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    if "CSP" in uservars and "CSP" in datavarlist:
        varlist.append(outvar(varname = "SnowNumberMix", longname = "Snow Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CSP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    if "CPP" in uservars and "CPP" in datavarlist:
        varlist.append(outvar(varname = "PrisNumberMix", longname = "Pristine Ice Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CPP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    if "CAP" in uservars and "CAP" in datavarlist:
        varlist.append(outvar(varname = "AggNumberMix", longname = "Aggregate Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CAP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    if "CGP" in uservars and "CGP" in datavarlist:
        varlist.append(outvar(varname = "GraupelNumberMix", longname = "Graupel Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CGP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    if "CHP" in uservars and "CHP" in datavarlist:
        varlist.append(outvar(varname = "HailNumberMix", longname = "Hail Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CHP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    if "CIPP" in uservars and "CIPP" in datavarlist:
        varlist.append(outvar(varname = "PlateNumberMix", longname = "Plate Ice Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CIPP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    if "CICP" in uservars and "CICP" in datavarlist:
        varlist.append(outvar(varname = "ColumnIceNumberMix", longname = "Column Ice Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CICP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    if "CIDP" in uservars and "CIDP" in datavarlist:
        varlist.append(outvar(varname = "DendriteNumberMix", longname = "Dendritic Snowflake Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CIDP", invar = None, unitfactor = 1, decnum = 3, units = "# kg**-1"))
    if "RCP" in uservars and "RCP" in datavarlist:
        varlist.append(outvar(varname = "CloudMassMix", longname = "Cloud Droplet Mixing Ratio", stdname = "cloud_liquid_water_mixing_ratio", vartype = "3d", nptype = float32, ramsname = "RCP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    if "RDP" in uservars and "RDP" in datavarlist:
        varlist.append(outvar(varname = "DrizzleMassMix", longname = "Drizzle Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RDP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    if "RRP" in uservars and "RRP" in datavarlist:
        varlist.append(outvar(varname = "RainMassMix", longname = "Rain Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RRP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    if "RSP" in uservars and "RSP" in datavarlist:
        varlist.append(outvar(varname = "SnowMassMix", longname = "Snow Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RSP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    if "RPP" in uservars and "RPP" in datavarlist:
        varlist.append(outvar(varname = "PrisMassMix", longname = "Pristine Ice Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RPP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    if "RAP" in uservars and "RAP" in datavarlist:
        varlist.append(outvar(varname = "AggMassMix", longname = "Aggregate Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RAP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    if "RGP" in uservars and "RGP" in datavarlist:
        varlist.append(outvar(varname = "GraupelMassMix", longname = "Graupel Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RGP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    if "RHP" in uservars and "RHP" in datavarlist:
        varlist.append(outvar(varname = "HailMassMix", longname = "Hail Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RHP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    if "RIPP" in uservars and "RIPP" in datavarlist:
        varlist.append(outvar(varname = "PlateMassMix", longname = "Plate Ice Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RIPP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    if "RICP" in uservars and "RICP" in datavarlist:
        varlist.append(outvar(varname = "ColumnIceMassMix", longname = "Column Ice Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RICP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))
    if "RICP" in uservars and "RICP" in datavarlist:
        varlist.append(outvar(varname = "DendriteMassMix", longname = "Dendritic Snowflake Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "RIDP", invar = None, unitfactor = "gperkg", decnum = 3, units = "g kg**-1"))

    #Hydrometeor Energies
    if "Q2" in uservars and "Q2" in datavarlist:
        varlist.append(outvar(varname = "RainEnergy", longname = "Internal Energy of Rain (J/kg Hydrometeor)", stdname = None, vartype = "3d", ramsname = "Q2", invar = None, unitfactor = 1, decnum = 3, units = "J kg**-1"))
    if "Q6" in uservars and "Q6" in datavarlist:
        varlist.append(outvar(varname = "GraupelEnergy", longname = "Internal Energy of Graupel (J/kg Hydrometeor)", stdname = None, vartype = "3d", ramsname = "Q6", invar = None, unitfactor = 1, decnum = 3, units = "J kg**-1"))
    if "Q7" in uservars and "Q7" in datavarlist:
        varlist.append(outvar(varname = "HailEnergy", longname = "Internal Energy of Hail (J/kg Hydrometeor)", stdname = None, vartype = "3d", ramsname = "Q7", invar = None, unitfactor = 1, decnum = 3, units = "J kg**-1"))

    #Hydrometeor growth budgets
    if "AGGREGATET" in uservars and "AGGREGATET" in datavarlist:
        varlist.append(outvar(varname = "Aggregation", longname = "Total Aggregation Between Analysis File Times", stdname = None, vartype = "3d", nptype = float32, ramsname = "AGGREGATET", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "AGGRPRISSNOWT" in uservars and "AGGRPRISSNOWT" in datavarlist:
        varlist.append(outvar(varname = "AggSnowPris", longname = "Total Aggregation from Collisions Between Pristine Ice and Snow Between Analysis File Times", stdname = None, vartype = "3d", nptype = float32, ramsname = "AGGRPRISSNOWT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "AGGRSELFSNOWT" in uservars and "AGGRSELFSNOWT" in datavarlist:
        varlist.append(outvar(varname = "AggSelfSnow", longname = "Total Aggregation from Collisions Between Aggregates and Snow Between Analysis File Times", stdname = None, vartype = "3d", nptype = float32, ramsname = "AGGRSELFSNOWT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "AGGRSELFPRIST" in uservars and "AGGRSELFSNOWT" in datavarlist:
        varlist.append(outvar(varname = "AggSelfPris", longname = "Total Aggregation from Collisions Between Aggregates and Pristine Ice Between Analysis File Times", stdname = None, vartype = "3d", nptype = float32, ramsname = "AGGRSELFPRIST", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "LATHEATFRZ" in uservars and "LATHEATFRZ" in datavarlist:
        varlist.append(outvar(varname = "LatentHeatFreezing", longname = "Rate of Temperature Change from Freezing/Melting", stdname = None, vartype = "3d", ramsname = "LATHEATFRZ", invar = None, unitfactor = "pertimestep", decnum = 6, units = "K s**-1"))
    if "LATHEATFRZT" in uservars and "LATHEATFRZT" in datavarlist:
        varlist.append(outvar(varname = "LatentHeatFreezingTotal", longname = "Total Temperature Change from Freezing/Melting Between Analysis File Times", stdname = None, vartype = "3d", ramsname = "LATHEATFRZT", invar = None, unitfactor = 1, decnum = 3, units = "K s**-1"))
    if "LATHEATVAP" in uservars and "LATHEATVAP" in datavarlist:
        varlist.append(outvar(varname = "LatentHeatCondensation", longname = "Rate of Temperature Change from Condensation/Evaporation and Deposition/Sublimation", stdname = None, vartype = "3d", ramsname = "LATHEATVAP", invar = None, unitfactor = "pertimestep", decnum = 6, units = "K s**-1"))
    if "LATHEATVAPT" in uservars and "LATHEATVAPT" in datavarlist:
        varlist.append(outvar(varname = "LatentHeatFreezingTotal", longname = "Total Temperature Change from Condensation/Evaporation and Deposition/Sublimation Between Analysis File Times", stdname = None, vartype = "3d", ramsname = "LATHEATVAPT", invar = None, unitfactor = 1, decnum = 3, units = "K s**-1"))
    if "CLD2RAINT" in uservars and "CLD2RAINT" in datavarlist:
        varlist.append(outvar(varname = "CollectedCloudWater", longname = "Rate of Collection of Cloud Droplets by Rain", stdname = None, vartype = "3d", ramsname = "CLD2RAINT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "ICE2RAINT" in uservars and "ICE2RAINT" in datavarlist:
        varlist.append(outvar(varname = "IceCollectionMelting", longname = "Rate of Ice Melting from Collection of Rain", stdname = None, vartype = "3d", ramsname = "ICE2RAINT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "RAIN2PRT" in uservars and "RAIN2PRT" in datavarlist:
        varlist.append(outvar(varname = "PrisRainFreezing", longname = "Rate of Ice Growth from Collisions of Rain and Pristine Ice", stdname = None, vartype = "3d", ramsname = "RAIN2PRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "RAIN2SNT" in uservars and "RAIN2SNT" in datavarlist:
        varlist.append(outvar(varname = "SnowRainRiming", longname = "Rate of Snow Growth from Riming of Rain", stdname = None, vartype = "3d", ramsname = "RAIN2SNT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    if "RAIN2AGT" in uservars and "RAIN2AGT" in datavarlist:
        varlist.append(outvar(varname = "AggRainRiming", longname = "Rate of Aggregate Growth from Riming of Rain", stdname = None, vartype = "3d", ramsname = "RAIN2AGT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "RAIN2GRT" in uservars and "RAIN2GRT" in datavarlist:
        varlist.append(outvar(varname = "GraupelRainRiming", longname = "Rate of Graupel Growth from Riming of Rain", stdname = None, vartype = "3d", ramsname = "RAIN2GRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "RAIN2HAT" in uservars and "RAIN2HAT" in datavarlist:
        varlist.append(outvar(varname = "HailRainRiming", longname = "Rate of Hail Growth from Riming of Rain", stdname = None, vartype = "3d", ramsname = "RAIN2HAT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "RIMECLDSNOWT" in uservars and "RIMECLDSNOWT" in datavarlist:
        varlist.append(outvar(varname = "SnowCloudRiming", longname = "Rate of Snow Growth from Riming of Cloud Droplets", stdname = None, vartype = "3d", ramsname = "RIMECLDSNOWT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "RIMECLDAGGRT" in uservars and "RIMECLDAGGRT" in datavarlist:
        varlist.append(outvar(varname = "AggCloudRiming", longname = "Rate of Aggregate Growth from Riming of Cloud Droplets", stdname = None, vartype = "3d", ramsname = "RIMECLDAGGRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "RIMECLDGRAUT" in uservars and "RIMECLDGRAUT" in datavarlist:
        varlist.append(outvar(varname = "GraupelCloudRiming", longname = "Rate of Graupel Growth from Riming of Cloud Droplets", stdname = None, vartype = "3d", ramsname = "RIMECLDGRAUT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "RIMECLDHAILT" in uservars and "RIMECLDHAILT" in datavarlist:
        varlist.append(outvar(varname = "HailCloudRiming", longname = "Rate of Hail Growth from Riming of Cloud Droplets", stdname = None, vartype = "3d", ramsname = "RIMECLDHAILT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "RIMECLDT" in uservars and "RIMECLDT" in datavarlist:
        varlist.append(outvar(varname = "TotalCloudRiming", longname = "Rate of Total Ice Growth from Riming of Cloud Droplets", stdname = None, vartype = "3d", ramsname = "RIMECLDT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "RAIN2ICET" in uservars and "RAIN2ICET" in datavarlist:
        varlist.append(outvar(varname = "TotalRainRiming", longname = "Rate of Total Ice Growth from Riming of Rain", stdname = None, vartype = "3d", ramsname = "RAIN2ICET", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "VAPCLDT" in uservars and "VAPCLDT" in datavarlist:
        varlist.append(outvar(varname = "CloudVapGrowth", longname = "Rate of Cloud Droplet Growth from Condensation", stdname = None, vartype = "3d", ramsname = "VAPCLDT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "VAPDRIZT" in uservars and "VAPDRIZT" in datavarlist:
        varlist.append(outvar(varname = "DrizzleVapGrowth", longname = "Rate of Drizzle Growth from Condensation", stdname = None, vartype = "3d", ramsname = "VAPDRIZT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "VAPRAINT" in uservars and "VAPRAINT" in datavarlist:
        varlist.append(outvar(varname = "RainVapGrowth", longname = "Rate of Rain Growth from Condensation", stdname = None, vartype = "3d", ramsname = "VAPRAINT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "VAPLIQT" in uservars and "VAPLIQT" in datavarlist:
        varlist.append(outvar(varname = "TotalLiquidVapGrowth", longname = "Rate of Total Liquid Growth from Condensation", stdname = None, vartype = "3d", ramsname = "VAPLIQT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "VAPPRIST" in uservars and "VAPPRIST" in datavarlist:
        varlist.append(outvar(varname = "PrisVapGrowth", longname = "Rate of Pristine Ice Growth from Vapor Deposition", stdname = None, vartype = "3d", ramsname = "VAPPRIST", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    if "VAPSNOWT" in uservars and "VAPSNOWT" in datavarlist:
        varlist.append(outvar(varname = "SnowVapGrowth", longname = "Rate of Snow Growth from Vapor Deposition", stdname = None, vartype = "3d", ramsname = "VAPSNOWT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    if "VAPAGGRT" in uservars and "VAPAGGRT" in datavarlist:
        varlist.append(outvar(varname = "AggVapGrowth", longname = "Rate of Aggregate Growth from Vapor Deposition", stdname = None, vartype = "3d", ramsname = "VAPAGGRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    if "VAPGRAUT" in uservars and "VAPGRAUT" in datavarlist:
        varlist.append(outvar(varname = "GraupelVapGrowth", longname = "Rate of Graupel Growth from Vapor Deposition", stdname = None, vartype = "3d", ramsname = "VAPGRAUT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    if "VAPHAILT" in uservars and "VAPHAILT" in datavarlist:
        varlist.append(outvar(varname = "HailVapGrowth", longname = "Rate of Hail Growth from Vapor Deposition", stdname = None, vartype = "3d", ramsname = "VAPHAILT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    if "VAPICET" in uservars and "VAPICET" in datavarlist:
        varlist.append(outvar(varname = "TotalIceVapGrowth", longname = "Rate of Total Ice Growth from Vapor Deposition", stdname = None, vartype = "3d", ramsname = "VAPICET", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**1"))
    if "EVAPCLDT" in uservars and "EVAPCLDT" in datavarlist:
        varlist.append(outvar(varname = "CloudEvap", longname = "Rate of Cloud Droplet Mass Loss from Evaporation", stdname = None, vartype = "3d", ramsname = "EVAPCLDT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "EVAPDRIZT" in uservars and "EVAPDRIZT" in datavarlist:
        varlist.append(outvar(varname = "DrizzleEvap", longname = "Rate of Drizzle Mass Loss from Evaporation", stdname = None, vartype = "3d", ramsname = "EVAPCLDT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "EVAPRAINT" in uservars and "EVAPRAINT" in datavarlist:
        varlist.append(outvar(varname = "RainEvap", longname = "Rate of Rain Mass Loss from Evaporation", stdname = None, vartype = "3d", ramsname = "EVAPRAINT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "EVAPLIQT" in uservars and "EVAPLIQT" in datavarlist:
        varlist.append(outvar(varname = "TotalLiquidEvap", longname = "Rate of Total Liquid Water Mass Loss from Evaporation", stdname = None, vartype = "3d", ramsname = "EVAPLIQT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "EVAPPRIST" in uservars and "EVAPPRIST" in datavarlist:
        varlist.append(outvar(varname = "PristineEvap", longname = "Rate of Pristine Ice Mass Loss from Sublimation", stdname = None, vartype = "3d", ramsname = "EVAPPRIST", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "EVAPSNOWT" in uservars and "EVAPSNOWT" in datavarlist:
        varlist.append(outvar(varname = "SnowEvap", longname = "Rate of Snow Mass Loss from Sublimation", stdname = None, vartype = "3d", ramsname = "EVAPPRIST", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "EVAPAGGRT" in uservars and "EVAPAGGRT" in datavarlist:
        varlist.append(outvar(varname = "AggEvap", longname = "Rate of Aggregate Mass Loss from Sublimation", stdname = None, vartype = "3d", ramsname = "EVAPAGGRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "EVAPGRAUT" in uservars and "EVAPGRAUT" in datavarlist:
        varlist.append(outvar(varname = "GraupelEvap", longname = "Rate of Graupel Mass Loss from Sublimation", stdname = None, vartype = "3d", ramsname = "EVAPGRAUT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "EVAPHAILT" in uservars and "EVAPHAILT" in datavarlist:
        varlist.append(outvar(varname = "HailEvap", longname = "Rate of Hail Mass Loss from Sublimation", stdname = None, vartype = "3d", ramsname = "EVAPHAILT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "MELTPRIST" in uservars and "MELTPRIST" in datavarlist:
        varlist.append(outvar(varname = "PrisMelting", longname = "Rate of Pristine Ice Melting", stdname = None, vartype = "3d", ramsname = "MELTPRIST", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "MELTSNOWT" in uservars and "MELTSNOWT" in datavarlist:
        varlist.append(outvar(varname = "SnowMelting", longname = "Rate of Snow Melting", stdname = None, vartype = "3d", ramsname = "MELTSNOWT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "MELTAGGRT" in uservars and "MELTAGGRT" in datavarlist:
        varlist.append(outvar(varname = "AggMelting", longname = "Rate of Aggregate Melting", stdname = None, vartype = "3d", ramsname = "MELTAGGRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "MELTGRAUT" in uservars and "MELTGRAUT" in datavarlist:
        varlist.append(outvar(varname = "GraupelMelting", longname = "Rate of Graupel Melting", stdname = None, vartype = "3d", ramsname = "MELTGRAUT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "MELTHAILT" in uservars and "MELTHAILT" in datavarlist:
        varlist.append(outvar(varname = "HailMelting", longname = "Rate of Hail Melting", stdname = None, vartype = "3d", ramsname = "MELTHAILT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "MELTICET" in uservars and "MELTICET" in datavarlist:
        varlist.append(outvar(varname = "TotalIceMelting", longname = "Rate of Total Ice Melting", stdname = None, vartype = "3d", ramsname = "MELTICET", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    
    #Freezing Processes
    if "INUCCONTRT" in uservars and "INUCCONTRT" in datavarlist:
        varlist.append(outvar(varname = "ContactFreezing", longname = "Rate of Contact Freezing", stdname = None, vartype = "3d", ramsname = "INUCCONTRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "INUCHOMRT" in uservars and "INUCHOMRT" in datavarlist:
        varlist.append(outvar(varname = "HomogeneousFreezing", longname = "Rate of Homogeneous Freezing", stdname = None, vartype = "3d", ramsname = "INUCHOMRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "INUCIFNRT" in uservars and "INUCIFNRT" in datavarlist:
        varlist.append(outvar(varname = "ImmersionFreezing", longname = "Rate of Condensation/Immersion Freezing", stdname = None, vartype = "3d", ramsname = "INUCIFNRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "INUCHAZRT" in uservars and "INUCHAZRT" in datavarlist:
        varlist.append(outvar(varname = "HazeFreezing", longname = "Rate of Haze Droplet Freezing", stdname = None, vartype = "3d", ramsname = "INUCHAZRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    
    #Dust Nucleation
    if "DUST1CLDRT" in uservars and "DUST1CLDRT" in datavarlist:
        varlist.append(outvar(varname = "Dust1CloudNucleation", longname = "Rate of Nucleation of Cloud Droplets from Dust Mode 1", stdname = None, vartype = "3d", ramsname = "DUST1CLDRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "DUST2CLDRT" in uservars and "DUST2CLDRT" in datavarlist:
        varlist.append(outvar(varname = "Dust2CloudNucleation", longname = "Rate of Nucleation of Cloud Droplets from Dust Mode 2", stdname = None, vartype = "3d", ramsname = "DUST2CLDRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "DUST1DRZRT" in uservars and "DUST1DRZRT" in datavarlist:
        varlist.append(outvar(varname = "Dust1DrizzleNucleation", longname = "Rate of Nucleation of Drizzle from Dust Mode 1", stdname = None, vartype = "3d", ramsname = "DUST1DRZRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    if "DUST2DRZRT" in uservars and "DUST2DRZRT" in datavarlist:
        varlist.append(outvar(varname = "Dust2DrizzleNucleation", longname = "Rate of Nucleation of Drizzle from Dust Mode 2", stdname = None, vartype = "3d", ramsname = "DUST2CLDRT", invar = None, unitfactor = "mgperkgpersec", decnum = 3, units = "mg kg**-1 s**-1"))
    
    #Aerosol Concentrations and Processes
    if "CCCNP" in uservars and "CCCNP" in datavarlist:
        varlist.append(outvar(varname = "CCNNumberMix", longname = "Cloud Condensation Nuclei Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CCCNP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "CIFNP" in uservars and "CIFNP" in datavarlist:
        varlist.append(outvar(varname = "INPNumberMix", longname = "Ice Nucleating Particle Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CIFNP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "GCCNP" in uservars and "GCCNP" in datavarlist:
        varlist.append(outvar(varname = "GCCNNumberMix", longname = "Giant CCN Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "GCCNP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "CN1NP" in uservars and "CN1NP" in datavarlist:
        varlist.append(outvar(varname = "CN1NumberMix", longname = "CCN Mode 1 Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CN1NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "CN2NP" in uservars and "CN2NP" in datavarlist:
        varlist.append(outvar(varname = "CN2NumberMix", longname = "CCN Mode 2 Number Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CN2NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "MD1NP" in uservars and "MD1NP" in datavarlist:
        varlist.append(outvar(varname = "SubDustNumberMix", longname = "Sub-Micron Dust Number Concentration", stdname = None, vartype = "3d", ramsname = "MD1NP", invar = None, unitfactor = 10**-6, decnum = 5, units = "# mg**-1"))
    if "MD2NP" in uservars and "MD2NP" in datavarlist:
        varlist.append(outvar(varname = "SuperDustNumberMix", longname = "Super-Micron Dust Number Concentration", stdname = None, vartype = "3d", ramsname = "MD2NP", invar = None, unitfactor = 10**-6, decnum = 5, units = "# mg**-1"))
    if "ABC1NP" in uservars and "ABC1NP" in datavarlist:
        varlist.append(outvar(varname = "ABC1NumberMix", longname = "Absorbing Carbon Type 1 Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "ABC1NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "ABC2NP" in uservars and "ABC2NP" in datavarlist:
        varlist.append(outvar(varname = "ABC2NumberMix", longname = "Absorbing Carbon Type 2 Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "ABC2NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "SALT_FILM_NP" in uservars and "SALT_FILM_NP" in datavarlist:
        varlist.append(outvar(varname = "SaltFilmNumberMix", longname = "Sea Salt Film Drop Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "SALT_FILM_NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "SALT_JET_NP" in uservars and "SALT_JET_NP" in datavarlist:
        varlist.append(outvar(varname = "SaltJetNumberMix", longname = "Sea Salt Jet Drop Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "SALT_JET_NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "SALT_SPUM_NP" in uservars and "SALT_SPUM_NP" in datavarlist:
        varlist.append(outvar(varname = "SaltSpumeNumberMix", longname = "Sea Salt Spume Drop Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "SALT_SPUM_NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "REGEN_AERO1_NP" in uservars and "REGEN_AERO1_NP" in datavarlist:
        varlist.append(outvar(varname = "SubAeroRegenNumberMix", longname = "Sub-Micron Regenerated Aerosol Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "REGEN_AERO1_NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "REGEN_AERO2_NP" in uservars and "REGEN_AERO2_NP" in datavarlist:
        varlist.append(outvar(varname = "SuperAeroRegenNumberMix", longname = "Super-Micron Regenerated Aerosol Number Mixing Ratio", stdname = None, vartype = "3d", ramsname = "REGEN_AERO2_NP", invar = None, unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "CCCMP" in uservars and "CCCMP" in datavarlist:
        varlist.append(outvar(varname = "CCNMassMix", longname = "Cloud Condensation Nuclei Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CCCMP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "GCCMP" in uservars and "GCCMP" in datavarlist:
        varlist.append(outvar(varname = "GCCNMassMix", longname = "Giant CCN Mixing Ratio", stdname = None, vartype = "3d", ramsname = "GCCMP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "CN1MP" in uservars and "CN1MP" in datavarlist:
        varlist.append(outvar(varname = "CN1MassMix", longname = "CCN Mode 1 Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CN1MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "CN2MP" in uservars and "CN2MP" in datavarlist:
        varlist.append(outvar(varname = "CN2MassMix", longname = "CCN Mode 2 Mixing Ratio", stdname = None, vartype = "3d", nptype = float32, ramsname = "CN2MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "MD1MP" in uservars and "MD1MP" in datavarlist:
        varlist.append(outvar(varname = "SubDustMassMix", longname = "Sub-Micron Dust Mixing Ratio", stdname = None, vartype = "3d", ramsname = "MD1MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "MD2MP" in uservars and "MD2MP" in datavarlist:
        varlist.append(outvar(varname = "SuperDustMassMix", longname = "Super-Micron Dust Mixing Ratio", stdname = None, vartype = "3d", ramsname = "MD2MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "ABC1MP" in uservars and "ABC1MP" in datavarlist:
        varlist.append(outvar(varname = "ABC1MassMix", longname = "Absorbing Carbon Type 1 Mixing Ratio", stdname = None, vartype = "3d", ramsname = "ABC1MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "ABC2MP" in uservars and "ABC2MP" in datavarlist:
        varlist.append(outvar(varname = "ABC2MassMix", longname = "Absorbing Carbon Type 2 Mixing Ratio", stdname = None, vartype = "3d", ramsname = "ABC2MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "SALT_FILM_MP" in uservars and "SALT_FILM_MP" in datavarlist:
        varlist.append(outvar(varname = "SaltFilmMassMix", longname = "Sea Salt Film Drop Mixing Ratio", stdname = None, vartype = "3d", ramsname = "SALT_FILM_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "SALT_JET_MP" in uservars and "SALT_JET_MP" in datavarlist:
        varlist.append(outvar(varname = "SaltJetMassMix", longname = "Sea Salt Jet Drop Mixing Ratio", stdname = None, vartype = "3d", ramsname = "SALT_JET_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "SALT_SPUM_MP" in uservars and "SALT_SPUM_MP" in datavarlist:
        varlist.append(outvar(varname = "SaltSpumeMassMix", longname = "Sea Salt Spume Drop Mixing Ratio", stdname = None, vartype = "3d", ramsname = "SALT_SPUM_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "REGEN_AERO1_MP" in uservars and "REGEN_AERO1_MP" in datavarlist:
        varlist.append(outvar(varname = "SubAeroRegenMassMix", longname = "Sub-Micron Regenerated Aerosol Mixing Ratio", stdname = None, vartype = "3d", ramsname = "REGEN_AERO1_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "REGEN_AERO2_MP" in uservars and "REGEN_AERO2_MP" in datavarlist:
        varlist.append(outvar(varname = "SuperAeroRegenMassMix", longname = "Super-Micron Regenerated Aerosol Mixing Ratio", stdname = None, vartype = "3d", ramsname = "REGEN_AERO2_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "RIFNP" in uservars and "RIFNP" in datavarlist:
        varlist.append(outvar(varname = "INPMassMix", longname = "Ice Nucleating Particle Mixing Ratio", stdname = None, vartype = "3d", ramsname = "RIFNP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))

    #3D Aerosol Tracking Variables
    if "IFNNUCP" in uservars and "IFNNUCP" in datavarlist:
        varlist.append(outvar(varname = "INPNucNumberMix", longname = "Number Mixing Ratio of INP Already Nucleated", stdname = None, vartype = "3d", ramsname = "IFNNUCP", unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "IMMERCP" in uservars and "IMMERCP" in datavarlist:
        varlist.append(outvar(varname = "CloudINPNumberMix", longname = "Number Mixing Ratio of INP Within Cloud Droplets", stdname = None, vartype = "3d", ramsname = "IMMERCP", unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "IMMERRP" in uservars and "IMMERRP" in datavarlist:
        varlist.append(outvar(varname = "RainINPNumberMix", longname = "Number Mixing Ratio of Within Rain", stdname = None, vartype = "3d", ramsname = "IMMERRP", unitfactor = 10**-6, decnum = 3, units = "# mg**-1"))
    if "CNMCP" in uservars and "CNMCP" in datavarlist:
        varlist.append(outvar(varname = "CloudAerosolMix", longname = "Mixing Ratio of Total Aerosol in Cloud Droplets", stdname = None, vartype = "3d", ramsname = "CNMCP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "CNMDP" in uservars and "CNMDP" in datavarlist:
        varlist.append(outvar(varname = "DrizzleAerosolMix", longname = "Mixing Ratio of Total Aerosol in Drizzle", stdname = None, vartype = "3d", ramsname = "CNMDP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "CNMRP" in uservars and "CNMRP" in datavarlist:
        varlist.append(outvar(varname = "RainAerosolMix", longname = "Mixing Ratio of Total Aerosol in Rain", stdname = None, vartype = "3d", ramsname = "CNMRP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "CNMPP" in uservars and "CNMPP" in datavarlist:
        varlist.append(outvar(varname = "PrisAerosolMix", longname = "Mixing Ratio of Total Aerosol in Pristine Ice", stdname = None, vartype = "3d", ramsname = "CNMPP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "CNMSP" in uservars and "CNMSP" in datavarlist:
        varlist.append(outvar(varname = "SnowAerosolMix", longname = "Mixing Ratio of Total Aerosol in Snow", stdname = None, vartype = "3d", ramsname = "CNMSP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "CNMAP" in uservars and "CNMAP" in datavarlist:
        varlist.append(outvar(varname = "AggAerosolMix", longname = "Mixing Ratio of Total Aerosol in Aggregates", stdname = None, vartype = "3d", ramsname = "CNMAP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "CNMGP" in uservars and "CNMGP" in datavarlist:
        varlist.append(outvar(varname = "GraupelAerosolMix", longname = "Mixing Ratio of Total Aerosol in Graupel", stdname = None, vartype = "3d", ramsname = "CNMGP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "CNMHP" in uservars and "CNMHP" in datavarlist:
        varlist.append(outvar(varname = "HailAerosolMix", longname = "Mixing Ratio of Total Aerosol in Hail", stdname = None, vartype = "3d", ramsname = "CNMHP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DNMCP" in uservars and "DNMCP" in datavarlist:
        varlist.append(outvar(varname = "CloudDustMix", longname = "Mixing Ratio of Dust in Cloud Droplets", stdname = None, vartype = "3d", ramsname = "DNMCP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DNMDP" in uservars and "DNMDP" in datavarlist:
        varlist.append(outvar(varname = "DrizzleDustMix", longname = "Mixing Ratio of Dust in Drizzle", stdname = None, vartype = "3d", ramsname = "DNMDP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DNMRP" in uservars and "DNMRP" in datavarlist:
        varlist.append(outvar(varname = "RainDustMix", longname = "Mixing Ratio of Dust in Rain", stdname = None, vartype = "3d", ramsname = "DNMRP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DNMPP" in uservars and "DNMPP" in datavarlist:
        varlist.append(outvar(varname = "PrisDustMix", longname = "Mixing Ratio of Dust in Pristine Ice", stdname = None, vartype = "3d", ramsname = "DNMPP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DNMSP" in uservars and "DNMSP" in datavarlist:
        varlist.append(outvar(varname = "SnowDustMix", longname = "Mixing Ratio of Dust in Snow", stdname = None, vartype = "3d", ramsname = "DNMSP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DNMAP" in uservars and "DNMAP" in datavarlist:
        varlist.append(outvar(varname = "AggDustMix", longname = "Mixing Ratio of Dust in Aggregates", stdname = None, vartype = "3d", ramsname = "DNMAP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DNMGP" in uservars and "DNMGP" in datavarlist:
        varlist.append(outvar(varname = "GraupelDustMix", longname = "Mixing Ratio of Dust in Graupel", stdname = None, vartype = "3d", ramsname = "DNMGP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DNMHP" in uservars and "DNMHP" in datavarlist:
        varlist.append(outvar(varname = "HailDustMix", longname = "Mixing Ratio of Dust in Hail", stdname = None, vartype = "3d", ramsname = "DNMHP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DINCP" in uservars and "DINCP" in datavarlist:
        varlist.append(outvar(varname = "CloudIceNucDust", longname = "Mixing Ratio of Dust in Cloud Droplets from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINCP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DINDP" in uservars and "DINDP" in datavarlist:
        varlist.append(outvar(varname = "DrizzleIceNucDust", longname = "Mixing Ratio of Dust in Drizzle from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINDP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DINRP" in uservars and "DINRP" in datavarlist:
        varlist.append(outvar(varname = "RainIceNucDust", longname = "Mixing Ratio of Dust in Rain from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINRP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DINPP" in uservars and "DINPP" in datavarlist:
        varlist.append(outvar(varname = "PrisIceNucDust", longname = "Mixing Ratio of Dust in Pristine Ice from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINPP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DINSP" in uservars and "DINSP" in datavarlist:
        varlist.append(outvar(varname = "SnowIceNucDust", longname = "Mixing Ratio of Dust in Snow from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINSP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DINAP" in uservars and "DINAP" in datavarlist:
        varlist.append(outvar(varname = "AggIceNucDust", longname = "Mixing Ratio of Dust in Aggregates from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINAP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DINGP" in uservars and "DINGP" in datavarlist:
        varlist.append(outvar(varname = "GraupelIceNucDust", longname = "Mixing Ratio of Dust in Graupel from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINGP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "DINHP" in uservars and "DINHP" in datavarlist:
        varlist.append(outvar(varname = "HailIceNucDust", longname = "Mixing Ratio of Dust in Hail from Ice Nucleation", stdname = None, vartype = "3d", ramsname = "DINHP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "SNMCP" in uservars and "SNMCP" in datavarlist:
        varlist.append(outvar(varname = "CloudSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Cloud Droplets", stdname = None, vartype = "3d", ramsname = "SNMCP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "SNMDP" in uservars and "SNMDP" in datavarlist:
        varlist.append(outvar(varname = "DrizzleSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Drizzle", stdname = None, vartype = "3d", ramsname = "SNMDP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "SNMRP" in uservars and "SNMRP" in datavarlist:
        varlist.append(outvar(varname = "RainSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Rain", stdname = None, vartype = "3d", ramsname = "SNMRP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "SNMPP" in uservars and "SNMPP" in datavarlist:
        varlist.append(outvar(varname = "PrisSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Pristine Ice", stdname = None, vartype = "3d", ramsname = "SNMPP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "SNMSP" in uservars and "SNMSP" in datavarlist:
        varlist.append(outvar(varname = "SnowSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Snow", stdname = None, vartype = "3d", ramsname = "SNMSP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "SNMAP" in uservars and "SNMAP" in datavarlist:
        varlist.append(outvar(varname = "AggSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Aggregates", stdname = None, vartype = "3d", ramsname = "SNMAP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "SNMGP" in uservars and "SNMGP" in datavarlist:
        varlist.append(outvar(varname = "GraupelSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Graupel", stdname = None, vartype = "3d", ramsname = "SNMGP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "SNMHP" in uservars and "SNMHP" in datavarlist:
        varlist.append(outvar(varname = "HailSolubleAeroMix", longname = "Mixing Ratio of Soluble Aerosol in Hail", stdname = None, vartype = "3d", ramsname = "SNMHP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "RESOL_AERO1_MP" in uservars and "RESOL_AERO1_MP" in datavarlist:
        varlist.append(outvar(varname = "RegenSolubleSubAeroMix", longname = "Mixing Ratio of Sub-Micron Regenerated Soluble Aerosol", stdname = None, vartype = "3d", ramsname = "RESOL_AERO1_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))
    if "RESOL_AERO2_MP" in uservars and "RESOL_AERO2_MP" in datavarlist:
        varlist.append(outvar(varname = "RegenSolubleSuperAeroMix", longname = "Mixing Ratio of Super-Micron Regenerated Soluble Aerosol", stdname = None, vartype = "3d", ramsname = "RESOL_AERO2_MP", invar = None, unitfactor = "ugperkg", decnum = 5, units = "ug kg**-1"))

    #3D Precipitation Variables
    if "PCPVD" in uservars and "PCPVD" in datavarlist:
        varlist.append(outvar(varname = "DrizzlePrecip3d", longname = "3D Rate of Drizzle Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVD", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPVR" in uservars and "PCPVR" in datavarlist:
        varlist.append(outvar(varname = "RainPrecip3d", longname = "3D Rate of Rain Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVR", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPVP" in uservars and "PCPVP" in datavarlist:
        varlist.append(outvar(varname = "PrisPrecip3d", longname = "3D Rate of Pristine Ice Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVP", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPVS" in uservars and "PCPVS" in datavarlist:
        varlist.append(outvar(varname = "SnowPrecip3d", longname = "3D Rate of Snow Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVS", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPVA" in uservars and "PCPVA" in datavarlist:
        varlist.append(outvar(varname = "AggPrecip3d", longname = "3D Rate of Aggregate Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVA", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPVG" in uservars and "PCPVG" in datavarlist:
        varlist.append(outvar(varname = "GraupelPrecip3d", longname = "3D Rate of Graupel Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVG", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPVH" in uservars and "PCPVH" in datavarlist:
        varlist.append(outvar(varname = "HailPrecip3d", longname = "3D Rate of Hail Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVH", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPVIP" in uservars and "PCPVIP" in datavarlist:
        varlist.append(outvar(varname = "PlatePrecip3d", longname = "3D Rate of Plate Ice Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVIP", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPVIC" in uservars and "PCPVIC" in datavarlist:
        varlist.append(outvar(varname = "ColumnIcePrecip3d", longname = "3D Rate of Column Ice Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVIC", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPVID" in uservars and "PCPVID" in datavarlist:
        varlist.append(outvar(varname = "DendritePrecip3d", longname = "3D Rate of Dendritic Snowflake Precipitation (mm/hr Liquid Equivalent)", stdname = None, vartype = "3d", ramsname = "PCPVIP", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))

    #Default Momentum Budget Variables
    if "WP_BUOY_THETA" in uservars and "WP_BUOY_THETA" in datavarlist:
        varlist.append(outvar(varname = "ThetaBuoyDefault", longname = "Default Model Vertical Acceleration from Potential Temperature Difference from Reference Potential Temperature", stdname = None, vartype = "3d", ramsname = "WP_BUOY_THETA", invar = None, unitfactor = "defaultmomentumbudget", decnum = 3, units = "m s**-2"))
    if "WP_BUOY_COND" in uservars and "WP_BUOY_COND" in datavarlist:
        varlist.append(outvar(varname = "CondBuoyDefault", longname = "Default Model Vertical Acceleration from Condensate Loading", stdname = None, vartype = "3d", ramsname = "WP_BUOY_COND", invar = None, unitfactor = "defaultmomentumbudget", decnum = 3, units = "m s**-2"))
    if "WP_ADVDIF" in uservars and "WP_ADVDIF" in datavarlist:
        varlist.append(outvar(varname = "AdvDifDefault", longname = "Default Model Vertical Acceleration from Advection and Diffusion", stdname = None, vartype = "3d", ramsname = "WP_ADVDIF", invar = None, unitfactor = "defaultmomentumbudget", decnum = 3, units = "m s**-2"))

    #Tracer Variables
    for i in range(1, 51):
        if f"TRACERP{str(i).zfill(3)}" in uservars and f"TRACERP{str(i).zfill(3)}" in datavarlist:
            varlist.append(outvar(varname = f"Tracer{str(i).zfill(3)}Mix", longname = f"Number Mixing Ratio of Tracer Number {str(i).zfill(3)}", stdname = None, nptype = float64, vartype = "3d", ramsname = f"TRACERP{str(i).zfill(3)}", unitfactor = 10**-6, decnum = 6, units = "# mg**-1"))

    #2D Grid Variables
    if "GLAT" in uservars and "GLAT" in datavarlist:
        varlist.append(outvar(varname = "Lat2D", longname = "2D Latitude Field from Grid", stdname = "grid_latitude", vartype = "2d", ramsname = "GLAT", invar = None, unitfactor = 1, decnum = 3, units = "degree"))
    if "GLON" in uservars and "GLON" in datavarlist:
        varlist.append(outvar(varname = "Lon2D", longname = "2D Longitude Field from Grid", stdname = "grid_longitude", vartype = "2d", ramsname = "GLON", invar = None, unitfactor = 1, decnum = 3, units = "degree"))
    if "TOPT" in uservars and "TOPT" in datavarlist:
        varlist.append(outvar(varname = "Topo", longname = "Topographic Height Above Mean Sea Level", stdname = "surface_altitude", vartype = "2d", ramsname = "TOPT", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    if "TOPZO" in uservars and "TOPZO" in datavarlist:
        varlist.append(outvar(varname = "TopoRoughness", longname = "Topographic Roughness Length", stdname = None, vartype = "2d", ramsname = "TOPZO", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    
    #2D Surface File Input Variables
    if "SOIL_TEXT" in uservars and "SOIL_TEXT" in datavarlist:
        varlist.append(outvar(varname = "SoilType", longname = "Soil Type of Current Soil Level", stdname = "soil_type", vartype = "2dSoil", ramsname = "SOIL_TEXT", invar = None, unitfactor = 1, decnum = 0, units = "Unitless"))
    if "LEAF_CLASS" in uservars and "LEAF_CLASS" in datavarlist:
        varlist.append(outvar(varname = "VegClass", longname = "Vegetation Type over Current Patch", stdname = "land_cover", vartype = "2dPatch", ramsname = "LEAF_CLASS", invar = None, unitfactor = 1, decnum = 0, units = "Unitless"))
    if "PATCH_AREA" in uservars and "PATCH_AREA" in datavarlist:
        varlist.append(outvar(varname = "PatchArea", longname = "Fraction of Grid Cell Occupied by Current Patch", stdname = None, vartype = "2dPatch", ramsname = "PATCH_AREA", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "VEG_NDVIP" in uservars and "VEG_NDVIP" in datavarlist:
        varlist.append(outvar(varname = "PastNDVI", longname = "Past Vegetation Normalized Difference Vegetation Index over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "VEG_NDVIP", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "VEG_NDVIC" in uservars and "VEG_NDVIC" in datavarlist:
        varlist.append(outvar(varname = "NDVI", longname = "Vegetation Normalized Difference Vegetation Index over Current Patch", stdname = "normalized_difference_vegetation_index", vartype = "2dPatch", ramsname = "VEG_NDVIC", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "VEG_NDVIF" in uservars and "VEG_NDVIF" in datavarlist:
        varlist.append(outvar(varname = "FutureNDVI", longname = "Future Vegetation Normalized Difference Vegetation Index over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "VEG_NDVIF", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    
    #2D Surface Characteristics 
    if "PATCH_ROUGH" in uservars and "PATCH_ROUGH" in datavarlist:
        varlist.append(outvar(varname = "PatchRough", longname = "Patch Roughness Length", stdname = None, vartype = "2dPatch", ramsname = "PATCH_ROUGH", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    if "SOIL_ROUGH" in uservars and "SOIL_ROUGH" in datavarlist:
        varlist.append(outvar(varname = "SoilRough", longname = "Soil Type Contribution to Roughness Length over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "SOIL_ROUGH", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    if "VEG_ROUGH" in uservars and "VEG_ROUGH" in datavarlist:
        varlist.append(outvar(varname = "VegRough", longname = "Vegetation Type Contribution to Roughness Length over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "VEG_ROUGH", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    if "VEG_FRACAREA" in uservars and "VEG_FRACAREA" in datavarlist:
        varlist.append(outvar(varname = "VegFrac", longname = "Fraction of Current Patch Covered by Vegetation", stdname = "vegetation_area_fraction", vartype = "2dPatch", ramsname = "VEG_FRACAREA", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    if "VEG_LAI" in uservars and "VEG_LAI" in datavarlist:
        varlist.append(outvar(varname = "VegLeafArea", longname = "Green Leaf Area Index (m^2 of Green Leaves to m^2 of Surface)", stdname = "leaf_area_index", vartype = "2dPatch", ramsname = "VEG_LAI", invar = None, unitfactor = 1, decnum = 3, units = "m**2 m**-2"))
    if "VEG_TAI" in uservars and "VEG_TAI" in datavarlist:
        varlist.append(outvar(varname = "VegTotalArea", longname = "Total Vegetation Area Index (m^2 of Green Leaves+Dead Leaves+Stems+ Trunks to m^2 of Surface)", stdname = None, vartype = "2dPatch", ramsname = "VEG_TAI", invar = None, unitfactor = 1, decnum = 3, units = "m**2 m**-2"))
    if "VEG_HEIGHT" in uservars and "VEG_HEIGHT" in datavarlist:
        varlist.append(outvar(varname = "VegHeight", longname = "Vegetation Height Above Ground", stdname = None, vartype = "2dPatch", ramsname = "VEG_HEIGHT", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    if "VEG_ALBEDO" in uservars and "VEG_ALBEDO" in datavarlist:
        varlist.append(outvar(varname = "VegAlbedo", longname = "Vegetation Broadband Albedo", stdname = None, vartype = "2dPatch", ramsname = "VEG_ALBEDO", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "SOIL_WATER" in uservars and "SOIL_WATER" in datavarlist:
        varlist.append(outvar(varname = "SoilMoisture", longname = "Volumetric Soil Moisture in Current Soil Level", stdname = "volume_fraction_of_condensed_water_in_soil", vartype = "2dSoil", ramsname = "SOIL_WATER", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "SOIL_ENERGY" in uservars and "SOIL_ENERGY" in datavarlist:
        varlist.append(outvar(varname = "SoilEnergy", longname = "Soil Energy With Respect to Freezing Within Current Soil Level", stdname = None, vartype = "2dSoil", ramsname = "SOIL_ENERGY", invar = None, unitfactor = 1, decnum = 3, units = "J kg**-1"))
    if "SFCWATER_NLEV" in uservars and "SFCWATER_NLEV" in datavarlist:
        varlist.append(outvar(varname = "SnowLevelCount", longname = "Number of Snow Levels over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "SFCWATER_NLEV", invar = None, unitfactor = 1, decnum = 0, units = "Unitless"))
    if "SFCWATER_MASS" in uservars and "SFCWATER_MASS" in datavarlist:
        varlist.append(outvar(varname = "SnowLevelMass", longname = "Snow Mass Within Current Snow Level", stdname = None, vartype = "2dSnow", ramsname = "SFCWATER_MASS", invar = None, unitfactor = 1, decnum = 3, units = "kg m**-2"))
    if "SFCWATER_ENERGY" in uservars and "SFCWATER_ENERGY" in datavarlist:
        varlist.append(outvar(varname = "SnowLevelEnergy", longname = "Snow Energy With Respect to Freezing Within Current Snow Level", stdname = None, vartype = "2dSnow", invar = None, ramsname = "SFCWATER_ENERGY", unitfactor = 1, decnum = 3, units = "J m**-3"))
    if "SFCWATER_DEPTH" in uservars and "SFCWATER_DEPTH" in datavarlist:
        varlist.append(outvar(varname = "SnowLevelDepth", longname = "Snow Depth Within Current Snow Level", stdname = None, vartype = "2dSnow", ramsname = "SFCWATER_DEPTH", invar = None, unitfactor = 1, decnum = 3, units = "m"))
    if "USTAR" in uservars and "USTAR" in datavarlist:
        varlist.append(outvar(varname = "UStar", longname = "Friction Velocity for Surface Fluxes over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "USTAR", invar = None, unitfactor = 1, decnum = 3, units = "m s**-1"))
    if "TSTAR" in uservars and "TSTAR" in datavarlist:
        varlist.append(outvar(varname = "TStar", longname = "Scaling Temperature for Surface Fluxes over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "TSTAR", invar = None, unitfactor = 1, decnum = 3, units = "K"))
    if "RSTAR" in uservars and "RSTAR" in datavarlist:
        varlist.append(outvar(varname = "RStar", longname = "Scaling Vapor Mixing Ratio for Surface Fluxes over Current Patch", stdname = None, vartype = "2dPatch", ramsname = "RSTAR", invar = None, unitfactor = 1000, decnum = 3, units = "g kg**-1"))
    if "STOM_RESIST" in uservars and "STOM_RESIST" in datavarlist:
        varlist.append(outvar(varname = "StomatalResistance", longname = "Vegetation Stomatal Resistance", stdname = None, vartype = "2dPatch", ramsname = "STOM_RESIST", invar = None, unitfactor = 1, decnum = 2, units = "s m**-1"))
    if "VEG_WATER" in uservars and "VEG_WATER" in datavarlist:
        varlist.append(outvar(varname = "VegSurfaceWater", longname = "Liquid Water Depth on Vegetation Surface", stdname = "canopy_water_amount", vartype = "2dPatch", ramsname = "VEG_WATER", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    if "VEG_TEMP" in uservars and "VEG_TEMP" in datavarlist:
        varlist.append(outvar(varname = "VegTemp", longname = "Vegetation Temperature", stdname = "canopy_temperature", vartype = "2dPatch", ramsname = "VEG_TEMP", invar = None, unitfactor = 1, decnum = 3, units = "K"))
    if "CAN_RVAP" in uservars and "CAN_RVAP" in datavarlist:
        varlist.append(outvar(varname = "CanopyVaporMix", longname = "Vapor Mixing Ratio in Vegetation Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "CAN_RVAP", invar = None, unitfactor = 1000, decnum = 3, units = "g kg**-1"))
    if "CAN_TEMP" in uservars and "CAN_TEMP" in datavarlist:
        varlist.append(outvar(varname = "CanopyTemp", longname = "Temperature in Vegetation Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "CAN_TEMP", invar = None, unitfactor = 1, decnum = 3, units = "K"))
    if "GROUND_RSAT" in uservars and "GROUND_RSAT" in datavarlist:
        varlist.append(outvar(varname = "GroundSatVaporMix", longname = "Saturation Vapor Mixing Ratio of Soil/Snow Surface (Used for Snow Sublimation and Dew Formation)", stdname = None, vartype = "2dPatch", ramsname = "GROUND_RSAT", invar = None, unitfactor = 1000, decnum = 3, units = "g kg**-1"))
    if "GROUND_RVAP" in uservars and "GROUND_RVAP" in datavarlist:
        varlist.append(outvar(varname = "SoilVaporMix", longname = "Vapor Mixing Ratio of Soil Surface (Used for Evaporation from Soil When No Snow Present)", stdname = None, vartype = "2dPatch", ramsname = "GROUND_RVAP", invar = None, unitfactor = 1000, decnum = 3, units = "g kg**-1"))

    #2D SiB Variables
    if "RCO2P" in uservars and "RCO2P" in datavarlist:
        varlist.append(outvar(varname = "CO2Mix", longname = "CO2 Mass Fraction (Divide by 1.5172*10**-6 to get CO2 ppm)", stdname = "mass_fraction_of_carbon_dioxide_in_air", vartype = "2dPatch", ramsname = "RCO2P", invar = None, unitfactor = 1, decnum = 7, units = "kg kg**-1"))
    if "SNOW1" in uservars and "SNOW1" in datavarlist:
        varlist.append(outvar(varname = "VegSnow", longname = "Mass of Snow Covering Vegetation", stdname = "canopy_snow_amount", vartype = "2dPatch", ramsname = "SNOW1", invar = None, unitfactor = 1, decnum = 4, units = "kg m**-2"))
    if "SNOW2" in uservars and "SNOW2" in datavarlist:
        varlist.append(outvar(varname = "SoilSnow", longname = "Mass of Snow Covering Bare Ground", stdname = "surface_snow_amount", vartype = "2dPatch", ramsname = "SNOW2", invar = None, unitfactor = 1, decnum = 4, units = "kg m**-2"))
    if "CAPAC1" in uservars and "CAPAC1" in datavarlist:
        varlist.append(outvar(varname = "VegStore", longname = "Vegetation Liquid Store", stdname = None, vartype = "2dPatch", ramsname = "CAPAC1", invar = None, unitfactor = 1, decnum = 4, units = "kg m**-2"))
    if "CAPAC2" in uservars and "CAPAC2" in datavarlist:
        varlist.append(outvar(varname = "SoilStore", longname = "Soil Surface Liquid Store", stdname = None, vartype = "2dPatch", ramsname = "CAPAC2", invar = None, unitfactor = 1, decnum = 4, units = "kg m**-2"))
    if "PCO2AP" in uservars and "PCO2AP" in datavarlist:
        varlist.append(outvar(varname = "CanopyCO2", longname = "Canopy Air Space CO2 Partial Pressure", stdname = None, vartype = "2dPatch", ramsname = "PCO2AP", invar = None, unitfactor = 1, decnum = 2, units = "Pa"))
    if "CO2FLX" in uservars and "CO2FLX" in datavarlist:
        varlist.append(outvar(varname = "CO2Flux", longname = "CO2 Flux from Canopy Air Space to Atmosphere", stdname = None, vartype = "2dPatch", ramsname = "CO2FLX", invar = None, unitfactor = 1, decnum = 4, units = "mol m**-2 s**-1"))
    if "SFCSWA" in uservars and "SFCSWA" in datavarlist:
        varlist.append(outvar(varname = "SIBSurfaceAlbedo", longname = "Broadband Surface Albedo of Current Patch", stdname = None, vartype = "2dPatch", ramsname = "SFCSWA", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "UPLWRF" in uservars and "UPLWR" in datavarlist:
        varlist.append(outvar(varname = "SIBLongwaveUpsrf", longname = "Upwelling Longwave Radiation From the Surface over the Current Patch", stdname = None, vartype = "2dPatch", ramsname = "UPLWRF", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "ASSIMN" in uservars and "ASSIMN" in datavarlist:
        varlist.append(outvar(varname = "CO2Uptake", longname = "Uptake of CO2 by Canopy Plants", stdname = None, vartype = "2dPatch", ramsname = "ASSIMN", invar = None, unitfactor = 1, decnum = 3, units = "umol m**-2 s**-1"))
    if "RESPG" in uservars and "RESPG" in datavarlist:
        varlist.append(outvar(varname = "GroundRespiration", longname = "Ground Respiration Flux", stdname = None, vartype = "2dPatch", ramsname = "RESPGF", invar = None, unitfactor = 1, decnum = 3, units = "umol m**-2 s**-1"))
    if "RSTFAC1" in uservars and "RSTFAC1" in datavarlist:
        varlist.append(outvar(varname = "LeafResistStress", longname = "Leaf Surface Humidity Resistance Stress", stdname = None, vartype = "2dPatch", ramsname = "RSTFAC1", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "RSTFAC2" in uservars and "RSTFAC2" in datavarlist:
        varlist.append(outvar(varname = "SoilResistStress", longname = "Soil Moisture Resistance Stress", stdname = None, vartype = "2dPatch", ramsname = "RSTFAC2", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "RSTFAC3" in uservars and "RSTFAC3" in datavarlist:
        varlist.append(outvar(varname = "TempResistStress", longname = "Temperature Resistance Stress", stdname = None, vartype = "2dPatch", ramsname = "RSTFAC3", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "ECT" in uservars and "ECT" in datavarlist:
        varlist.append(outvar(varname = "Transpiration", longname = "Latent Heat Flux into Canopy Air Space from Transpiration", stdname = "transpiration_flux", vartype = "2dPatch", ramsname = "ECT", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "ECI" in uservars and "ECI" in datavarlist:
        varlist.append(outvar(varname = "CanopyIntercept", longname = "Latent Heat Flux into Canopy Air Space from Evaporation of Water on Vegetation Surface", stdname = "tendency_of_canopy_water_amount_from_evaporation_of_intercepted_precipitation", vartype = "2dPatch", ramsname = "ECI", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "EGI" in uservars and "EGI" in datavarlist:
        varlist.append(outvar(varname = "GroundIntercept", longname = "Latent Heat Flux into Canopy Air Space from Evaporation of Water on Soil Surface", stdname = None, vartype = "2dPatch", ramsname = "EGI", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "EGS" in uservars and "EGS" in datavarlist:
        varlist.append(outvar(varname = "GroundEvap", longname = "Latent Heat Flux into Canopy Air Space from Evaporation of Soil Moisture", stdname = "water_evaporation_flux_from_soil", vartype = "2dPatch", ramsname = "EGS", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "HC" in uservars and "HC" in datavarlist:
        varlist.append(outvar(varname = "VegSensHeatFlux", longname = "Sensible Heat Flux into Canopy Air Space from Vegetation", stdname = None, vartype = "2dPatch", ramsname = "HC", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "HG" in uservars and "HG" in datavarlist:
        varlist.append(outvar(varname = "SoilSensHeatFlux", longname = "Sensible Heat Flux into Canopy Air Space from Soil Surface", stdname = None, vartype = "2dPatch", ramsname = "HG", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "RA" in uservars and "RA" in datavarlist:
        varlist.append(outvar(varname = "CanopyToAtmosResistance", longname = "Aerodynamic Resistance from Canopy Air Space to Atmosphere", stdname = None, vartype = "2dPatch", ramsname = "RA", invar = None, unitfactor = 1, decnum = 2, units = "s m**-1"))
    if "RB" in uservars and "RB" in datavarlist:
        varlist.append(outvar(varname = "VegToCanopyResistance", longname = "Aerodynamic Resistance from Vegetation Surface to Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "RB", invar = None, unitfactor = 1, decnum = 2, units = "s m**-1"))
    if "RC" in uservars and "RC" in datavarlist:
        varlist.append(outvar(varname = "TotalCanopyResistance", longname = "Total Aerodynamic Resistance of Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "RC", invar = None, unitfactor = 1, decnum = 2, units = "s m**-1"))
    if "RD" in uservars and "RD" in datavarlist:
        varlist.append(outvar(varname = "SoilToCanopyResistance", longname = "Aerodynamic Resistance from Soil Surface to Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "RD", invar = None, unitfactor = 1, decnum = 2, units = "s m**-1"))
    if "ROFF" in uservars and "ROFF" in datavarlist:
        varlist.append(outvar(varname = "Runoff", longname = "Depth of Surface and Subsurface Runoff", stdname = None, vartype = "2dPatch", ramsname = "ROFF", invar = None, unitfactor = 1, decnum = 2, units = "mm"))
    if "GREEN" in uservars and "GREEN" in datavarlist:
        varlist.append(outvar(varname = "GreenFrac", longname = "Greenness Fraction of Vegetation (Area of Green Leaves to Total Plant Area)", stdname = None, vartype = "2dPatch", ramsname = "GREEN", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "APAR" in uservars and "APAR" in datavarlist:
        varlist.append(outvar(varname = "PhotoAbsorb", longname = "Fraction of Photosynthetically Active Radiation Absorbed by Vegetation", stdname = None, vartype = "2dPatch", ramsname = "APAR", invar = None, unitfactor = 1, decnum = 4, units = "Unitless"))
    if "VENTMF" in uservars and "VENTMF" in datavarlist:
        varlist.append(outvar(varname = "VentilationMassFlux", longname = "Ventilation Mass Flux", stdname = None, vartype = "2dPatch", ramsname = "VENTMF", invar = None, unitfactor = 1, decnum = 4, units = "kg m**-2 s**-1"))
    if "PCO2C" in uservars and "PCO2C" in datavarlist:
        varlist.append(outvar(varname = "ChloroCO2Pressure", longname = "Partial Pressure of CO2 in Chloroplasts", stdname = None, vartype = "2dPatch", ramsname = "PCO2C", invar = None, unitfactor = 1, decnum = 2, units = "Pa"))
    if "PCO2I" in uservars and "PCO2I" in datavarlist:
        varlist.append(outvar(varname = "LeafInternalCO2Pressure", longname = "Partial Pressure of CO2 in Internal Leaf Space", stdname = None, vartype = "2dPatch", ramsname = "PCO2I", invar = None, unitfactor = 1, decnum = 2, units = "Pa"))
    if "PCO2S" in uservars and "PCO2S" in datavarlist:
        varlist.append(outvar(varname = "LeafSurfaceCO2Pressure", longname = "Partial Pressure of CO2 at Leaf Surface", stdname = None, vartype = "2dPatch", ramsname = "PCO2S", invar = None, unitfactor = 1, decnum = 2, units = "Pa"))
    if "EA" in uservars and "EA" in datavarlist:
        varlist.append(outvar(varname = "CanopyVaporPres", longname = "Water Vapor Pressure in Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "EA", invar = None, unitfactor = 1, decnum = 3, units = "hPa"))
    if "EM" in uservars and "EM" in datavarlist:
        varlist.append(outvar(varname = "ReferenceVaporPres", longname = "Reference Level Water Vapor Pressure", stdname = None, vartype = "2dPatch", ramsname = "EM", invar = None, unitfactor = 1, decnum = 3, units = "hPa"))
    if "RHA" in uservars and "RHA" in datavarlist:
        varlist.append(outvar(varname = "CanopyHumidity", longname = "Relative Humidity of Canopy Air Space", stdname = None, vartype = "2dPatch", ramsname = "RHA", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "RADVBC" in uservars and "RADVBC" in datavarlist:
        varlist.append(outvar(varname = "VisibleDiffuse", longname = "Visible Diffuse Radiation", stdname = "surface_diffuse_downwelling_shortwave_flux_in_air", vartype = "2dPatch", ramsname = "RADVBC", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "RADVDC" in uservars and "RADVDC" in datavarlist:
        varlist.append(outvar(varname = "VisibleDirect", longname = "Visible Direct Radiation", stdname = "surface_direct_downwelling_shortwave_flux_in_air", vartype = "2dPatch", ramsname = "RADVDC", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "RADNBC" in uservars and "RADNBC" in datavarlist:
        varlist.append(outvar(varname = "NIRDiffuse", longname = "Near-Infrared Diffuse Radiation", stdname = None, vartype = "2dPatch", ramsname = "RADNBC", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "RADNDC" in uservars and "RADNDC" in datavarlist:
        varlist.append(outvar(varname = "NIRDirect", longname = "Near-Infrared Direct Radiation", stdname = None, vartype = "2dPatch", ramsname = "RADNDC", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "PSY" in uservars and "PSY" in datavarlist:
        varlist.append(outvar(varname = "PsyConstant", longname = "Psychrometric Constant", stdname = None, vartype = "2dPatch", ramsname = "PSY", invar = None, unitfactor = 1, decnum = 3, units = "hPa deg**-1"))

    #2D Precipitation Variables
    if "PCPRD" in uservars and "PCPRD" in datavarlist:
        varlist.append(outvar(varname = "DrizzlePrecipRate", longname = "Rate of Drizzle Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRD", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPRR" in uservars and "PCPRR" in datavarlist:
        varlist.append(outvar(varname = "RainPrecipRate", longname = "Rate of Rain Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = "rainfall_rate", vartype = "2d", ramsname = "PCPRR", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPRP" in uservars and "PCPRP" in datavarlist:
        varlist.append(outvar(varname = "PrisPrecipRate", longname = "Rate of Pristine Ice Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRP", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPRS" in uservars and "PCPRS" in datavarlist:
        varlist.append(outvar(varname = "SnowPrecipRate", longname = "Rate of Snow Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRS", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPRA" in uservars and "PCPRA" in datavarlist:
        varlist.append(outvar(varname = "AggPrecipRate", longname = "Rate of Aggregate Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRA", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPRG" in uservars and "PCPRG" in datavarlist:
       varlist.append(outvar(varname = "GraupelPrecipRate", longname = "Rate of Graupel Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRG", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPRH" in uservars and "PCPRH" in datavarlist:
        varlist.append(outvar(varname = "HailPrecipRate", longname = "Rate of Hail Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRH", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPRIP" in uservars and "PCPRIP" in datavarlist:
        varlist.append(outvar(varname = "PlatePrecipRate", longname = "Rate of Plate Ice Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRIP", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPRIC" in uservars and "PCPRIC" in datavarlist:
        varlist.append(outvar(varname = "ColumnIcePrecipRate", longname = "Rate of Column Ice Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRIC", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "PCPRID" in uservars and "PCPRID" in datavarlist:
        varlist.append(outvar(varname = "DendritePrecipRate", longname = "Rate of Dendritic Snowflake Precipitation at the Surface (mm/hr Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "PCPRID", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "ACCPD" in uservars and "ACCPD" in datavarlist:
        varlist.append(outvar(varname = "DrizzlePrecipTotal", longname = "Total Accumulated Drizzle Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPD", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    if "ACCPR" in uservars and "ACCPR" in datavarlist:
        varlist.append(outvar(varname = "RainPrecipTotal", longname = "Total Accumulated Rain Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPR", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    if "ACCPP" in uservars and "ACCPP" in datavarlist:
        varlist.append(outvar(varname = "PrisPrecipTotal", longname = "Total Accumulated Pristine Ice Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPP", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    if "ACCPS" in uservars and "ACCPS" in datavarlist:
        varlist.append(outvar(varname = "SnowPrecipTotal", longname = "Total Accumulated Snow Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPS", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    if "ACCPA" in uservars and "ACCPA" in datavarlist:
        varlist.append(outvar(varname = "AggPrecipTotal", longname = "Total Accumulated Aggregate Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPA", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    if "ACCPG" in uservars and "ACCPG" in datavarlist:
        varlist.append(outvar(varname = "GraupelPrecipTotal", longname = "Total Accumulated Graupel Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPG", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    if "ACCPH" in uservars and "ACCPH" in datavarlist:
        varlist.append(outvar(varname = "HailPrecipTotal", longname = "Total Accumulated Hail Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPH", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    if "ACCPIP" in uservars and "ACCPIP" in datavarlist:
        varlist.append(outvar(varname = "PlatePrecipTotal", longname = "Total Accumulated Plate Ice Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPIP", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    if "ACCPIC" in uservars and "ACCPIC" in datavarlist:
        varlist.append(outvar(varname = "ColumnIcePrecipTotal", longname = "Total Accumulated Column Ice Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPIC", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    if "ACCPID" in uservars and "ACCPID" in datavarlist:
        varlist.append(outvar(varname = "DendritePrecipTotal", longname = "Total Accumulated Dendritic Snowflake Precipitation at the Surface (mm Liquid Equivalent)", stdname = None, vartype = "2d", ramsname = "ACCPID", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    if "PCPG" in uservars and "PCPG" in datavarlist:
        varlist.append(outvar(varname = "LEAFPrecipRate", longname = "Instantaneous Mass Flux of Microphysics Precipitation (used by LEAF/SIB)", stdname = None, vartype = "2d", ramsname = "PCPG", invar = None, unitfactor = "pertimestep", decnum = 3, units = "kg m**-2"))
    if "QPCPG" in uservars and "QPCPG" in datavarlist:
        varlist.append(outvar(varname = "LEAFPrecipEnergy", longname = "Instantaneous Energy Flux of Microphysics Precipitation (used by LEAF/SIB)", stdname = None, vartype = "2d", ramsname = "QPCPG", invar = None, unitfactor = "pertimestep", decnum = 3, units = "W m**-2"))
    if "DPCPG" in uservars and "DPCPG" in datavarlist:
        varlist.append(outvar(varname = "LEAFPrecipDepth", longname = "Instantaneous Rate of Precipitation Depth (used by LEAF/SIB)", stdname = None, vartype = "2d", ramsname = "DPCPG", invar = None, unitfactor = "pertimestep", decnum = 3, units = "m s**-1"))
    
    #Aerosol Deposition
    if "ACCPDUST" in uservars and "ACCPDUST" in datavarlist:
        varlist.append(outvar(varname = "TotalDustDeposition", longname = "Total Deposition of Dust on the Surface", stdname = None, vartype = "2d", ramsname = "ACCPDUST", invar = None, unitfactor = 1, decnum = 3, units = "kg m**-2"))
    if "ACCPAERO" in uservars and "ACCPAERO" in datavarlist:
        varlist.append(outvar(varname = "TotalAerosolDeposition", longname = "Total Deposition of Aerosol on the Surface", stdname = None, vartype = "2d", ramsname = "ACCPAERO", invar = None, unitfactor = 1, decnum = 3, units = "kg m**-2"))
    if "PCPRDUST" in uservars and "PCPRDUST" in datavarlist:
        varlist.append(outvar(varname = "DustDepositionRate", longname = "Rate of Deposition of Dust on the Surface", stdname = None, vartype = "2d", ramsname = "PCPRDUST", invar = None, unitfactor = 1, decnum = 3, units = "kg m**-2 s**-1"))
    if "PCPRAERO" in uservars and "PCPRAERO" in datavarlist:
        varlist.append(outvar(varname = "AeroDepositionRate", longname = "Rate of Deposition of Aerosol on the Surface", stdname = None, vartype = "2d", ramsname = "PCPRAERO", invar = None, unitfactor = 1, decnum = 3, units = "kg m**-2"))

    #2D Radiation Variables
    if "RSHORT" in uservars and "RSHORT" in datavarlist:
        varlist.append(outvar(varname = "ShortwaveDownSrf", longname = "Surface Downwelling Shortwave Radiation", stdname = "surface_downwelling_shortwave_flux", vartype = "2d", ramsname = "RSHORT", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "RLONG" in uservars and "RLONG" in datavarlist:
        varlist.append(outvar(varname = "LongwaveDownSrf", longname = "Surface Downwelling Longwave Radiation", stdname = "surface_downwelling_longwave_flux", vartype = "2d", ramsname = "RLONG", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "RLONGUP" in uservars and "RLONGUP" in datavarlist:
        varlist.append(outvar(varname = "LongwaveUpSrf", longname = "Surface Upwelling Longwave Radiation", stdname = "surface_upwelling_longwave_flux", vartype = "2d", ramsname = "RLONGUP", invar = None, unitfactor = 1, decnum = 3, units = "W m**-2"))
    if "AODT" in uservars and "AODT" in datavarlist:
        varlist.append(outvar(varname = "AerosolOpticalDepth", longname = "Aerosol Optical Depth in Visible Light", stdname = "atmosphere_optical_thickness_due_to_ambient_aerosol_particles", vartype = "2d", ramsname = "AODT", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "ALBEDT" in uservars and "ALBEDT" in datavarlist:
        varlist.append(outvar(varname = "SurfaceAlbedo", longname = "Surface Broadband Albedo", stdname = "surface_albedo", vartype = "2d", ramsname = "ALBEDT", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    if "COSZ" in uservars and "COSZ" in datavarlist:
        varlist.append(outvar(varname = "SolarZenithAngle", longname = "Solar Zenith Angle", stdname = "surface_albedo", vartype = "2d", ramsname = "COSZ", invar = None, unitfactor = 1, decnum = 3, units = "Unitless"))
    
    #2D Surface Fluxes
    if "SFLUX_T" in uservars and "SFLUX_T" in datavarlist:
        varlist.append(outvar(varname = "SensibleHeatFlux", longname = "Surface Sensible Heat Flux", stdname = "surface_upward_sensible_heat_flux", vartype = "2d", ramsname = "SFLUX_T", invar = None, unitfactor = 1004, decnum = 3, units = "W m**-2"))
    if "SFLUX_R" in uservars and "SFLUX_R" in datavarlist:
        varlist.append(outvar(varname = "LatentHeatFlux", longname = "Surface Latent Heat Flux", stdname = "surface_upward_latent_heat_flux", vartype = "2d", ramsname = "SFLUX_R", invar = None, unitfactor = 2.5*10**6, decnum = 3, units = "W m**-2"))
    if "SFLUX_U" in uservars and "SFLUX_U" in datavarlist:
        varlist.append(outvar(varname = "UMomentumFlux", longname = "Surface U-Momentum Flux", stdname = None, vartype = "2d", ramsname = "SFLUX_U", invar = None, unitfactor = 1, decnum = 3, units = "Pa"))
    if "SFLUX_V" in uservars and "SFLUX_U" in datavarlist:
        varlist.append(outvar(varname = "VMomentumFlux", longname = "Surface V-Momentum Flux", stdname = None, vartype = "2d", ramsname = "SFLUX_V", invar = None, unitfactor = 1, decnum = 3, units = "Pa"))
    if "SFLUX_W" in uservars and "SFLUX_W" in datavarlist:
        varlist.append(outvar(varname = "WMomentumFlux", longname = "Surface W-Momentum Flux", stdname = None, vartype = "2d", ramsname = "SFLUX_W", invar = None, unitfactor = 1, decnum = 3, units = "Pa"))
    
    #2D Convective Parameterization Variables
    if "ACONPR" in uservars and "ACONPR" in datavarlist:
        varlist.append(outvar(varname = "ConvectivePrecipTotal", longname = "Total Accumulated Precipitation at the Surface from Convective Parameterization", stdname = "convective_precipitation_amount", vartype = "2d", ramsname = "ACONPR", invar = None, unitfactor = 1, decnum = 3, units = "mm"))
    if "CONPRR" in uservars and "CONPRR" in datavarlist:
        varlist.append(outvar(varname = "ConvectivePrecipRate", longname = "Rate of Precipitation at the Surface from Convective Parameterization", stdname = "convective_precipitation_rate", vartype = "2d", ramsname = "CONPRR", invar = None, unitfactor = "preciprate", decnum = 3, units = "mm hr**-1"))
    if "NCA" in uservars and "NCA" in datavarlist:
        varlist.append(outvar(varname = "ConvectiveTimesteps", longname = "Number of Timesteps at Horizontal Gridpoint Maintaining Convective Behavior", stdname = None, vartype = "2d", ramsname = "NCA", invar = None, unitfactor = 1, decnum = 3, units = 1))
    if "CONVGO" in uservars and "CONVGO" in datavarlist:
        varlist.append(outvar(varname = "ConvectiveInitCheck", longname = "Check if Pre-Convection Checks at Horizontal Gridpoint Satisfied", stdname = None, vartype = "2d", ramsname = "CONVGO", invar = None, unitfactor = 1, decnum = 3, units = 1))
    
    
    
    return varlist

# def datavarinit(varentry, ftime, coords):
#     pass

def gen_vardict(varlist):
    vdict = {}
    for entry in varlist:
        vdict[entry.ramsname] = entry
    return vdict