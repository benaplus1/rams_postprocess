import pickle
from varinit import outvar

def derivvarinit():
    dvardict = {} #This dictionary contains *all* possible derived variables. Based on conditions in derivedvars_2.py, only some will be calculated
    dvardict["TS"] = outvar(varname = "SrfTemp", longname = "Terrain-Following Air Temperature at the Lowest Model Level", stdname = None, offline = True, vartype = "2d", ramsname = "TS", invar = ["THETA", "PI"], unitfactor = 1, decnum = 3, units = "K")

    dvardict["PS"] = outvar(varname = "SrfPres", longname = "Terrain-Following Pressure at the Lowest Model Level", stdname = None, offline = True, vartype = "2d", ramsname = "PS", invar = ["PI"], unitfactor = 1, decnum = 3, units = "hPa")

    dvardict["MSLP"] = outvar(varname = "MSLP", longname = "Mean Sea Level Pressure", stdname = None, offline = True, vartype = "2d", ramsname = "MSLP", invar = ["PI", "THETA", "TOPT"], unitfactor = 1, decnum = 3, units = "hPa") #Note the 0.0065 here is a 6.5 C/km lapse rate, which is also used by REVU's MSLP calculation. It may be worth revisiting this with a different lapse rate formulation in future.
    
    dvardict["MSLT"] = outvar(varname = "MSLT", longname = "Mean Sea Level Temperature", stdname = None, offline = True, vartype = "2d", ramsname = "MSLT", invar = ["THETA", "TOPT"], unitfactor = 1, decnum = 3, units = "K") #Again, using a 6.5 C/km lapse rate to approximate sea level temperature. Need a better approach than a constant lapse rate, though. Testing in my lake-effect sims shows that 9.8 C/km actually works much better, so it might depend on the situation
    
    dvardict["RTC"] = outvar(varname = "CondMix", longname = "Total Condensate Mixing Ratio", stdname = None, offline = True, vartype = "3d", ramsname = "RTC", invar = ["RTP", "RV"], unitfactor = "gperkg", decnum = 3, units = "g kg**-1")
    
    dvardict["RTI"] = outvar(varname = "IceMassMix", longname = "Total Ice Mixing Ratio", stdname = None, offline = True, vartype = "3d", ramsname = "RTI", invar = ["RPP", "RSP", "RAP", "RGP", "RHP"], unitfactor = "gperkg", decnum = 3, units = "g kg**-1")
    
    dvardict["RTL"] = outvar(varname = "LiqMassMix", longname = "Total Liquid Mixing Ratio", stdname = None, offline = True, vartype = "3d", ramsname = "RTL", invar = ["RCP", "RDP", "RRP"], unitfactor = "gperkg", decnum = 3, units = "g kg**-1")
    
    dvardict["THV"] = outvar(varname = "ThetaV", longname = "Virtual Potential Temperature", stdname = None, offline = True, vartype = "3d", ramsname = "THV", invar = ["RV", "THETA"], unitfactor = 1, decnum = 3, units = "K")
    
    dvardict["THR"] = outvar(varname = "ThetaRho", longname = "Density Potential Temperature", stdname = None, offline = True, vartype = "3d", ramsname = "THR", invar = ["RV", "THETA", "RTP"], unitfactor = 1, decnum = 3, units = "K")
    
    dvardict["P"] = outvar(varname = "Pressure", longname = "Air Pressure", stdname = "air_pressure", offline = True, vartype = "3d", ramsname = "P", invar = ["PI"], unitfactor = 1, decnum = 3, units = "hPa")
    
    dvardict["T"] = outvar(varname = "Temperature", longname = "Air Temperature", stdname = "air_temperature", offline = True, vartype = "3d", ramsname = "T", invar = ["THETA", "PI"], unitfactor = 1, decnum = 3, units = "K")
    
    dvardict["TV"] = outvar(varname = "TempV", longname = "Virtual Temperature", stdname = "virtual_temperature", offline = True, vartype = "3d", ramsname = "TV", invar = ["THETA", "PI", "RV"], unitfactor = 1, decnum = 3, units = "K")
   
    dvardict["TR"] = outvar(varname = "TempRho", longname = "Density Temperature", stdname = None, offline = True, vartype = "3d", ramsname = "TR", invar = ["THETA", "PI", "RV", "RTP"], unitfactor = 1, decnum = 3, units = "K")
    
    dvardict["RHO"] = outvar(varname = "Rho", longname = "Air Density", stdname = "air_density", offline = True, vartype = "3d", ramsname ="RHO", invar = ["THETA", "PI", "RV"], unitfactor = 1, decnum = 4, units = "kg m**-3")

    dvardict["DIV"] = outvar(varname = "HorizDiv", longname = "Horizontal Divergence", stdname = None, offline = True, vartype = "3d", ramsname = "DIV", invar = ["UC", "VC"], unitfactor = 1, decnum = 7, units = "s**-1")
    
    dvardict["E"] = outvar(varname = "VaporPressure", longname = "Vapor Pressure", stdname = "water_vapor_partial_pressure_in_air", offline = True, vartype = "3d", ramsname = "E", invar = ["PI", "RV"], unitfactor = 1, decnum = 5, units = "hPa")
    
    dvardict["PWAT"] = outvar(varname = "PWAT", longname = "Precipitable Water (Total Column Water Vapor)", stdname = None, offline = True, vartype = "2d", ramsname = "PWAT", invar = ["THETA", "PI", "RV"], unitfactor = 1, decnum = 3, units = "mm")
    
    dvardict["TD"] = outvar(varname = "Dewpoint", longname = "Dew Point Temperature", stdname = "dew_point_temperature", offline = True, vartype = "3d", ramsname = "TD", invar = ["RV", "PI"], unitfactor = 1, decnum = 3, units = "K")
    
    dvardict["ZCT"] = outvar(varname = "CloudTopHeight", longname = "Cloud Top Height (Including Ice Clouds)", stdname = "cloud_top_altitude", offline = True, vartype = "2d", ramsname = "ZCT", invar = ["RCP", "RPP"], unitfactor = 1, decnum = 0, units = "m")
    
    dvardict["ZCB"] = outvar(varname = "CloudBaseHeight", longname = "Cloud Base Height (Including Ice Clouds)", stdname = "cloud_base_altitude", offline = True, vartype = "2d", ramsname = "ZCB", invar = ["RCP", "RPP"], unitfactor = 1, decnum = 0, units = "m")
    
    dvardict["PCT"] = outvar(varname = "CloudTopPressure", longname = "Pressure at Cloud Top", stdname = "air_pressure_at_cloud_top", offline = True, vartype = "2d", ramsname = "PCT", invar = ["RCP", "RPP", "PI"], unitfactor = 1, decnum = 2, units = "hPa")
    
    dvardict["PCB"] = outvar(varname = "CloudBasePressure", longname = "Pressure at Cloud Base", stdname = "air_pressure_at_cloud_base", offline = True, vartype = "2d", ramsname = "PCB", invar = ["RCP", "RPP", "PI"], unitfactor = 1, decnum = 2, units = "hPa")
    
    dvardict["WADVH"] = outvar(varname = "WAdv_H", longname = "Horizontal Advection of Vertical Velocity", stdname = None, offline = True, vartype = "3d", ramsname = "WADVH", invar = ["UC", "VC", "WC"], unitfactor = 1, decnum = 5, units = "m s**-2")
    
    dvardict["WADVV"] = outvar(varname = "WAdv_V", longname = "Vertical Advection of Vertical Velocity", stdname = None, offline = True, vartype = "3d", ramsname = "WADVV", invar = ["WC"], unitfactor = 1, decnum = 5, units = "m s**-2")
    
    dvardict["VOX"] = outvar(varname = "XVort", longname = "Relative Vorticity in the +X (Eastward) Direction", stdname = "atmosphere_x_relative_vorticity", offline = True, vartype = "3d", ramsname = "VOX", invar = ["VC", "WC"], unitfactor = 1, decnum = 7, units = "s**-1")

    dvardict["VOY"] = outvar(varname = "YVort", longname = "Relative Vorticity in the +Y (Northward) Direction", stdname = "atmosphere_y_relative_vorticity", offline = True, vartype = "3d", ramsname = "VOY", invar = ["UC", "WC"], unitfactor = 1, decnum = 7, units = "s**-1")

    dvardict["VOZR"] = outvar(varname = "ZVortRel", longname = "Relative Vorticity in the +Z (Upward) Direction", stdname = "atmosphere_upward_relative_vorticity", offline = True, vartype = "3d", ramsname = "VOZR", invar = ["UC", "VC"], unitfactor = 1, decnum = 7, units = "s**-1")

    dvardict["VOZA"] = outvar(varname = "ZVortAbs", longname = "Absolute Vorticity in the +Z (Upward) Direction", stdname = "atmosphere_upward_absolute_vorticity", offline = True, vartype = "3d", ramsname = "VOZA", invar = ["UC", "VC"], unitfactor = 1, decnum = 7, units = "s**-1")

    dvardict["SL"] = outvar(varname = "LiqSatFrac", longname = "Water Vapor Pressure Fraction of Equilibrium Vapor Pressure with Respect to Liquid", stdname = None, offline = True, vartype = "3d", ramsname = "SL", invar = ["THETA", "PI", "RV"], unitfactor = 1, decnum = 4, units = "Unitless")

    dvardict["DVL"] = outvar(varname = "LiqVaporDeficit", longname = "Difference Between Equilibrium Vapor Pressure With Respect to Liquid and Current Water Vapor Pressure", stdname = "water_vapor_saturation_deficit_in_air", offline = True, vartype = "3d", ramsname = "DVL", invar = ["THETA", "PI", "RV"], unitfactor = 1, decnum = 3, units = "hPa")

    dvardict["SI"] = outvar(varname = "IceSatFrac", longname = "Water Vapor Pressure Fraction of Equilibrium Vapor Pressure with Respect to Ice", stdname = None, offline = True, vartype = "3d", ramsname = "SI", invar = ["THETA", "PI", "RV"], unitfactor = 1, decnum = 4, units = "Unitless")

    dvardict["DVI"] = outvar(varname = "IceVaporDeficit", longname = "Difference Between Equilibrium Vapor Pressure with Respect to Ice and Current Water Vapor Pressure", stdname = None, offline = True, vartype = "3d", ramsname = "DVI", invar = ["THETA", "PI", "RV"], unitfactor = 1, decnum = 3, units = "hPa")

    dvardict["RVSL"] = outvar(varname = "SuperLiqMix", longname = "Mixing Ratio Difference from Equilibrium Water Vapor Mixing Ratio with Respect to Liquid", stdname = None, offline = True, vartype = "3d", ramsname = "RVSL", invar = ["THETA", "PI", "RV"], unitfactor = "gperkg", decnum = 6, units = "g kg**-1")

    dvardict["RVSI"] = outvar(varname = "SuperIceMix", longname = "Mixing Ratio Difference from Equilibrium Water Vapor Mixing Ratio With Respect to Ice", stdname = None, offline = True, vartype = "3d", ramsname = "RVSI", invar = ["THETA", "PI", "RV"], unitfactor = "gperkg", decnum = 6, units = "g kg**-1")

    #Hydrometeor Number Concentrations
    dvardict["NCP"] = outvar(varname = "CloudNum", longname = "Cloud Droplet Number Concentration", stdname = "number_concentration_of_cloud_liquid_water_particles_in_air", offline = True, vartype = "3d", ramsname = "NCP", invar = ["CCP", "DN0"], unitfactor = 1, decnum = 0, units = "m**-3")

    dvardict["NDP"] = outvar(varname = "DrizzleNum", longname = "Drizzle Drop Number Concentration", stdname = None, offline = True, vartype = "3d", ramsname = "NDP", invar = ["CDP", "DN0"], unitfactor = 1, decnum = 0, units = "m**-3")
    
    dvardict["NRP"] = outvar(varname = "RainNum", longname = "Rain Drop Number Concentration", stdname = None, offline = True, vartype = "3d", ramsname = "NRP", invar = ["CRP", "DN0"], unitfactor = 1, decnum = 0, units = "m**-3")
    
    dvardict["NPP"] = outvar(varname = "PrisNum", longname = "Pristine Ice Number Concentration", stdname = None, offline = True, vartype = "3d", ramsname = "NPP", invar = ["CPP", "DN0"], unitfactor = 1, decnum = 0, units = "m**-3")

    dvardict["NSP"] = outvar(varname = "SnowNum", longname = "Snow Number Concentration", stdname = None, offline = True, vartype = "3d", ramsname = "NSP", invar = ["CSP", "DN0"], unitfactor = 1, decnum = 0, units = "m**-3")

    dvardict["NAP"] = outvar(varname = "AggNum", longname = "Aggregate Number Concentration", stdname = None, offline = True, vartype = "3d", ramsname = "NAP", invar = ["CAP", "DN0"], unitfactor = 1, decnum = 0, units = "m**-3")

    dvardict["NGP"] = outvar(varname = "GraupelNum", longname = "Graupel Number Concentration", stdname = None, offline = True, vartype = "3d", ramsname = "NGP", invar = ["CGP", "DN0"], unitfactor = 1, decnum = 0, units = "m**-3")

    dvardict["NHP"] = outvar(varname = "HailNum", longname = "Hail Number Concentration", stdname = None, offline = True, vartype = "3d", ramsname = "NHP", invar = ["CHP", "DN0"], unitfactor = 1, decnum = 0, units = "m**-3")

    dvardict["NIPP"] = outvar(varname = "PlateNum", longname = "Plate Ice Number Concentration", stdname = None, offline = True, vartype = "3d", ramsname = "NIPP", invar = ["CIPP", "DN0"], unitfactor = 1, decnum = 0, units = "m**-3")

    dvardict["NICP"] = outvar(varname = "ColumnNum", longname = "Column Ice Number Concentration", stdname = None, offline = True, vartype = "3d", ramsname = "NICP", invar = ["CICP", "DN0"], unitfactor = 1, decnum = 0, units = "m**-3")

    dvardict["NIDP"] = outvar(varname = "DendriteNum", longname = "Dendritic Snowflake Number Concentration", stdname = None, offline = True, vartype = "3d", ramsname = "NIDP", invar = ["CIDP", "DN0"], unitfactor = 1, decnum = 0, units = "m**-3")

    #Hydrometeor Water Contents
    dvardict["CWC"] = outvar(varname = "CloudWaterContent", longname = "Cloud Droplet Water Content", stdname = "mass_concentration_of_cloud_liquid_water_in_air", offline = True, vartype = "3d", ramsname = "CWC", invar = ["RCP", "DN0"], unitfactor = 1, decnum = 7, units = "kg m**-3")

    dvardict["DWC"] = outvar(varname = "DrizzleWaterContent", longname = "Drizzle Drop Water Content", stdname = None, offline = True, vartype = "3d", ramsname = "DWC", invar = ["RDP", "DN0"], unitfactor = 1, decnum = 7, units = "kg m**-3")

    dvardict["RWC"] = outvar(varname = "RainWaterContent", longname = "Rain Water Content", stdname = None, offline = True, vartype = "3d", ramsname = "RWC", invar = ["RRP", "DN0"], unitfactor = 1, decnum = 7, units = "kg m**-3")

    dvardict["PWC"] = outvar(varname = "PrisWaterContent", longname = "Pristine Ice Water Content", stdname = None, offline = True, vartype = "3d", ramsname = "PWC", invar = ["RPP", "DN0"], unitfactor = 1, decnum = 7, units = "kg m**-3")

    dvardict["SWC"] = outvar(varname = "SnowWaterContent", longname = "Snow Water Content", stdname = None, offline = True, vartype = "3d", ramsname = "SWC", invar = ["RSP", "DN0"], unitfactor = 1, decnum = 7, units = "kg m**-3")

    dvardict["AWC"] = outvar(varname = "AggWaterContent", longname = "Aggregate Water Content", stdname = None, offline = True, vartype = "3d", ramsname = "AWC", invar = ["RAP", "DN0"], unitfactor = 1, decnum = 7, units = "kg m**-3")

    dvardict["GWC"] = outvar(varname = "GraupelWaterContent", longname = "Graupel Water Content", stdname = None, offline = True, vartype = "3d", ramsname = "GWC", invar = ["RGP", "DN0"], unitfactor = 1, decnum = 7, units = "kg m**-3")

    dvardict["HWC"] = outvar(varname = "HailWaterContent", longname = "Hail Water Content", stdname = None, offline = True, vartype = "3d", ramsname = "HWC", invar = ["RHP", "DN0"], unitfactor = 1, decnum = 7, units = "kg m**-3")

    dvardict["IPWC"] = outvar(varname = "PlateWaterContent", longname = "Plate Ice Water Content", stdname = None, offline = True, vartype = "3d", ramsname = "IPWC", invar = ["RIPP", "DN0"], unitfactor = 1, decnum = 7, units = "kg m**-3")

    dvardict["ICWC"] = outvar(varname = "ColumnWaterContent", longname = "Column Ice Water Content", stdname = None, offline = True, vartype = "3d", ramsname = "ICWC", invar = ["RICP", "DN0"], unitfactor = 1, decnum = 7, units = "kg m**-3")

    dvardict["IDWC"] = outvar(varname = "DendriteWaterContent", longname = "Dendritic Snowflake Water Content", stdname = None, offline = True, vartype = "3d", ramsname = "IDWC", invar = ["RIDP", "DN0"], unitfactor = 1, decnum = 7, units = "kg m**-3")
    
    #Vertically Integrated Hydrometeor Amounts
    dvardict["VIC"] = outvar(varname = "VertIntCloud", longname = "Liquid Equivalent of Vertically Integrated Cloud Water", stdname = None, offline = True, vartype = "2d", ramsname = "VIC", invar = ["RCP", "DN0"], unitfactor = 1, decnum = 3, units = "mm")
        #Again, using the full sigma-z grid to get all cloud water up to the model top, not just the top of the user-chosen analysis domain

    dvardict["VID"] = outvar(varname = "VertIntDrizzle", longname = "Liquid Equivalent of Vertically Integrated Drizzle", stdname = None, offline = True, vartype = "2d", ramsname = "VID", invar = ["RDP", "DN0"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["VIR"] = outvar(varname = "VertIntRain", longname = "Liquid Equivalent of Vertically Integrated Rain", stdname = None, offline = True, vartype = "2d", ramsname = "VIR", invar = ["RRP", "DN0"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["VIP"] = outvar(varname = "VertIntPris", longname = "Liquid Equivalent of Vertically Integrated Pristine Ice", stdname = None, offline = True, vartype = "2d", ramsname = "VIP", invar = ["RPP", "DN0"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["VIS"] = outvar(varname = "VertIntSnow", longname = "Liquid Equivalent of Vertically Integrated Snow", stdname = None, offline = True, vartype = "2d", ramsname = "VIS", invar = ["RSP", "DN0"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["VIA"] = outvar(varname = "VertIntAgg", longname = "Liquid Equivalent of Vertically Integrated Aggregates", stdname = None, offline = True, vartype = "2d", ramsname = "VIA", invar = ["RAP", "DN0"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["VIG"] = outvar(varname = "VertIntGraupel", longname = "Liquid Equivalent of Vertically Integrated Graupel", stdname = None, offline = True, vartype = "2d", ramsname = "VIG", invar = ["RGP", "DN0"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["VIH"] = outvar(varname = "VertIntHail", longname = "Liquid Equivalent of Vertically Integrated Hail", stdname = None, offline = True, vartype = "2d", ramsname = "VIH", invar = ["RHP", "DN0"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["VIIP"] = outvar(varname = "VertIntPlate", longname = "Liquid Equivalent of Vertically Integrated Plate Ice", stdname = None, offline = True, vartype = "2d", ramsname = "VIIP", invar = ["RIPP", "DN0"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["VIIC"] = outvar(varname = "VertIntColumn", longname = "Liquid Equivalent of Vertically Integrated Column Ice", stdname = None, offline = True, vartype = "2d", ramsname = "VIIC", invar = ["RICP", "DN0"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["VIID"] = outvar(varname = "VertIntDendrite", longname = "Liquid Equivalent of Vertically Integrated Dendritic Snowflakes", stdname = None, offline = True, vartype = "2d", ramsname = "VIID", invar = ["RIDP", "DN0"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["VIL"] = outvar(varname = "VertIntLiq", longname = "Vertically Integrated Liquid", stdname = None, offline = True, vartype = "2d", ramsname = "VIL", invar = ["RCP", "RDP", "RRP"], unitfactor = 1, decnum = 2, units = "mm")

    dvardict["VII"] = outvar(varname = "VertIntIce", longname = "Vertically Integrated Ice", stdname = None, offline = True, vartype = "2d", ramsname = "VII", invar = ["RPP", "RSP", "RAP", "RGP", "RHP"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["VIT"] = outvar(varname = "VertIntCond", longname = "Verticaly Integrated Total Condensate", stdname = None, offline = True, vartype = "2d", ramsname = "VIT", invar = ["RTP", "RV"], unitfactor = 1, decnum = 4, units = "mm")
    
    #Hydrometeor Diameters
    #I don't need to use the gamma function for this, thankfully!
    dvardict["DCP"] = outvar(varname = "CloudDiam", longname = "Mass Mean Cloud Droplet Diameter", stdname = None, offline = True, vartype = "3d", ramsname = "DCP", invar = ["RCP", "CCP"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["DDP"] = outvar(varname = "DrizzleDiam", longname = "Mass Mean Drizzle Drop Diameter", stdname = None, offline = True, vartype = "3d", ramsname = "DDP", invar = ["RDP", "CDP"], unitfactor = 1, decnum = 3, units = "mm")

    dvardict["DRP"] = outvar(varname = "RainDiam", longname = "Mass Mean Rain Drop Diameter", stdname = None, offline = True, vartype = "3d", ramsname = "DRP", invar = ["RRP", "CRP"], unitfactor = 1, decnum = 2, units = "mm")

    dvardict["DPP"] = outvar(varname = "PrisDiam", longname = "Mass Mean Pristine Ice Diameter", stdname = None, offline = True, vartype = "3d", ramsname = "DPP", invar = ["RPP", "CPP"], unitfactor = 1, decnum = 4, units = "mm")

    dvardict["DSP"] = outvar(varname = "SnowDiam", longname = "Mass Mean Snow Diameter", stdname = None, offline = True, vartype = "3d", invar = ["RSP", "CSP"], unitfactor = 1, decnum = 3, units = "mm")

    dvardict["DAP"] = outvar(varname = "AggDiam", longname = "Mass Mean Aggregate Diameter", stdname = None, offline = True, vartype = "3d", invar = ["RAP", "CAP"], unitfactor = 1, decnum = 2, units = "mm")

    dvardict["DGP"] = outvar(varname = "GraupelDiam", longname = "Mass Mean Graupel Diameter", stdname = None, offline = True, vartype = "3d", invar = ["RGP", "CHP"], unitfactor = 1, decnum = 2, units = "mm")

    dvardict["DHP"] = outvar(varname = "HailDiam", longname = "Mass Mean Hail Diameter", stdname = None, offline = True, vartype = "3d", invar = ["RHP", "CHP"], unitfactor = 1, decnum = 2, units = "mm")
    
    #Reference State and Momentum Budget Calculations. Currently, can only be done for cartesian interpolation, unfortunately
    dvardict["RHOP"] = outvar(varname = "RhoPrime", longname = "Density Perturbation (Difference of Air Density at Current Grid Cell from Local Mean Air Density)", stdname = None, offline = True, vartype = "3d", ramsname = "RHOP", invar = ["THETA", "PI", "RV", "RTP"], unitfactor = 1, decnum = 5, units = "kg m**-3")

    dvardict["BUOY_RHO"] = outvar(varname = "RhoBuoy", longname = "Buoyancy Calculated from Density Perturbation", stdname = None, offline = True, vartype = "3d", ramsname = "BUOY_RHO", invar = ["THETA", "PI", "RV", "RTP"], unitfactor = 1, decnum = 5, units = "m s**-2")

    dvardict["BUOY_THETA"] = outvar(varname = "ThetaRhoBuoy", longname = "Buoyancy Calculated from Density Potential Temperature Perturbation", stdname = None, offline = True, vartype = "3d", ramsname = "BUOY_THETA", invar = ["THETA", "PI", "RV", "RTP"], unitfactor = 1, decnum = 5, units = "m s**_2")

    dvardict["BUOY_TEMP"] = outvar(varname = "ThermalBuoy", longname = "Buoyancy Due to Temperature Difference from Environment", stdname = None, offline = True, vartype = "3d", ramsname = "BUOY_TEMP", invar = ["THETA"], unitfactor = 1, decnum = 5, units = "m s**-2")

    dvardict["BUOY_VAP"] = outvar(varname = "VaporBuoy", longname = "Buoyancy Due to Decrease in Density from Water Vapor", stdname = None, offline = True, vartype = "3d", ramsname = "BUOY_VAP", invar = ["RV"], unitfactor = 1, decnum = 5, units = "m s**-2")

    dvardict["BUOY_COND"] = outvar(varname = "CondBuoy", longname = "Negative Buoyancy from Condensate Loading", stdname = None, offline = True, vartype = "3d", ramsname = "BUOY_COND", invar = ["RTP", "RV"], unitfactor = 1, decnum = 5, units = "m s**-2")

    dvardict["BUOY_PPRIME"] = outvar(varname = "PPrimeBuoy", longname = "Buoyancy Due to Horizontal Pressure Perturbations", stdname = None, offline = True, vartype = "3d", ramsname = "BUOY_PPRIME", invar = ["PI"], unitfactor = 1, decnum = 5, units = "m s**-2")

    dvardict["VPPGF"] = outvar(varname = "VPPGF", longname = "Acceleration due to Vertical Perturbation Pressure Gradient", stdname = None, offline = True, vartype = "3d", ramsname = "VPPGF", invar = ["THETA", "PI", "RV", "RTP"], unitfactor = 1, decnum = 5, units = "m s**-2")
    
    return dvardict

#This dictionary links the "verbose names" and the "rams names" for each derived variable. We need this because in the backend, I use the rams names for all derived variables. 
if __name__ == "__main__":
    # dvarnamedict = {"SrfTemp": "TS", "SrfPres": "PS", "MSLP": "MSLP", "MSLT": "MSLT", "CondMix": "RTC",
    # "IceMix": "RTI", "LiqMix": "RTL", "ThetaV": "THV", "ThetaRho": "THR", "Pressure": "P",
    # "Temperature": "T", "TempV": "TV", "TempRho": "TR", "Rho": "RHO", "VaporPressure": "E",
    # "PWAT": "PWAT", "Dewpoint": "TD", "CloudTopHeight": "ZCT", "CloudBaseHeight": "ZCB",
    # "CloudTopPressure": "PCT", "CloudBasePressure": "PCB", "WAdv_H": "WADVH",
    # "WAdv_V": "WADVV", "XVort": "VOX", "YVort": "VOY", "ZVortRel": "VOZR",
    # "ZVortAbs": "VOZA", "LiqSatFrac": "SL", "IceSatFrac": "SI", "LiqVaporDeficit": "DVL",
    # "IceVaporDeficit": "DVI", "SuperLiqMix": "RVSI", "SuperIceMix": "RVSI", "CloudNum": "NCP",
    # "DrizNum": "NDP", "RainNum": "NRP", "PrisNum": "NPP", "SnowNum": "NSP", "AggNum": "NAP",
    # "GraupelNum": "NGP", "HailNum": "NHP", "PlateNum": "NIPP", "ColumnNum": "NICP",
    # "DendriteNum": "NIDP", "CloudWaterContent": "CWC", "DrizzleWaterContent": "DWC",
    # "RainWaterContent": "RWC", "PrisWaterContent": "PWC", "SnowWaterContent": "SWC",
    # "AggWaterContent": "AWC", "GraupelWaterContent": "GWC", "HailWaterContent": "HWC",
    # "PlateWaterContent": "IPWC", "ColumnWaterContent": "ICWC", "DendriteWaterContent": "IDWC",
    # "VertIntCloud": "VIC", "VertIntDrizzle": "VID", "VertIntRain": "VIR",
    # "VertIntPris": "VIP", "VertIntSnow": "VIS", "VertIntAgg": "VIA", "VertIntGraupel": "VIG",
    # "VertIntHail": "VIH", "VertIntPlate": "VIIP", "VertIntColumn": "VIIC",
    # "VertIntDendrite": "VIID", "VertIntLiq": "VIL", "VertIntIce": "VII",
    # "CloudDiam": "DCP", "DrizzleDiam": "DDP", "RainDiam": "DRP", "PrisDiam": "DPP",
    # "SnowDiam": "DSP", "AggDiam": "DAP", "GraupelDiam": "DGP", "HailDiam": "DHP",
    # "RhoPrime": "RHOP", "RhoBuoy": "BUOY_RHO", "ThetaRhoBuoy": "BUOY_THETA",
    # "ThermalBuoy": "BUOY_TEMP", "VapBuoy": "BUOY_VAP", "CondBuoy": "BUOY_COND", "PPrimeBuoy": "BUOY_PPRIME", "VPPGF": "VPPGF"}

    dvardict = derivvarinit()
    dvarnamedict = {dvardict[key].varname:dvardict[key].ramsname for key in dvardict.keys()}
    print(dvarnamedict)

    with open("dvardictfile", "wb") as wfile:
        pickle.dump(dvarnamedict, wfile)
    print("Done!")