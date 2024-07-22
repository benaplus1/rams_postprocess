#!/usr/bin/env python
# coding: utf-8

# In[2]:


import xarray as xr
import pickle
import h5netcdf
import numpy as np
import pandas as pd
from time import perf_counter
from datetime import datetime, timedelta
from pandas import date_range
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from interplib import get_gprops, gen_coords_cart, get_zlvs_sigma, get_zact, gen_coords_sigma, gen_coords_pressure, get_cartw, get_sigw, get_presw, get_cartu, get_sigu, get_presu, get_cartv, get_sigv, get_presv, get_cartscl_multiprocess, get_presscl_multiprocess, get_sigscl_multiprocess
from varinit import outvar, fillvarlist, gen_vardict
from derivedvars import get_derivedvars

import os
def read_userfile(fpath):
    with open(fpath, "r") as userfile:
        mlist = userfile.read().splitlines()
        folderpath = mlist[mlist.index("RAMS Analysis File folder path")+1]
        postpath = mlist[mlist.index("Folder in which to put post-processed data")+1]
        filetype = mlist[mlist.index("Do you want to process Analysis files or Lite files? Enter *L* for Lite, *A* for Analysis")+1]
        ramsinpath = mlist[mlist.index("Path to RAMSIN for this simulation")+1]
        # project = input("Enter the project (idealized, escape, or les): ")
        # case = input("Enter the case/simulation to use: ")
        interptype = mlist[mlist.index("Type of vertical coordinates to use: Cartesian coordinates, Pressure coordinates, or Sigma coordinates")+1]
        atop = mlist[mlist.index("If using cartesian or sigma interpolation, the maximum altitude (in meters) you'd like to include in post-processing. leave blank to analyze the whole grid")+1]
        userpreslvs = mlist[mlist.index("If using pressure interpolation, the pressure levels in hPa you'd like the post-processed data on as comma-separated values")+1]
        instime = mlist[mlist.index("Analysis start time in YYYY-mm-dd-HHMMSS format")+1]
        inetime = mlist[mlist.index("Analysis end time in YYYY-mm-dd-HHMMSS format")+1]
        numworkers = mlist[mlist.index("Number of cores to use for post-processing")+1]
        ngrid = mlist[mlist.index("Grid Number to post-process")+1] #For now, can only do one grid at a time; if you want to analyze a fine and coarse grid, will need to run this twice with different ngrid
        uservars = mlist[mlist.index("List of RAMS variables you'd like to process, as comma-separated entries. Put *all* to process all available variables")+1]
        rnameflag = mlist[mlist.index("Use the *RAMS* variable names, or *verbose* variable names in post-processed NetCDF files?")+1]
        derivvars = mlist[mlist.index("List of derived variables you'd like to output, as comma-separated entries (a full list of available derived variables is avilable in comments at the top of derivedvars.py). Leave blank to not output any derived quantities. Enter *all* to output all derived quantities available for your output file. Enter *nomomentum* to output all variables except momentum budgets, which are quite slow to calculate")+1]
        ywininput = mlist[mlist.index("If outputting momentum budgets, the number of Y grid points used for horizontal averaging")+1]
        xwininput = mlist[mlist.index("If outputting momentum budgets, the number of X grid points used for horizontal averaging")+1]
        kernname = mlist[mlist.index("If outputting momentum budgets, the type of convolution kernel used for horizontal averaging (documentation is available in derivedvars.py)")+1]
        userdict = {"folderpath": folderpath, "postpath": postpath, "filetype": filetype, "ramsinpath": ramsinpath, "interptype": interptype, "atop": atop, "userpreslvs": userpreslvs, "instime": instime, "inetime": inetime, "numworkers": numworkers, "ngrid": ngrid, "uservars": uservars, "rnameflag": rnameflag, "derivvars": derivvars, "ywininput": ywininput, "xwininput": xwininput, "kernname": kernname}
    return userdict

def get_timedict(ramsinpath, ngrid, ftype):
    # Get grid properties from RAMS header file
    with open(ramsinpath) as f:
        mlist = f.read().splitlines()
        dtlong = float([st for st in mlist if st.lstrip().startswith("DTLONG")][0].lstrip("DTLONG = ").split(",")[0])
        filetime_ramsin_varname = "FRQSTATE" if ftype == "A" else "FRQLITE"
        afiletime = int([st for st in mlist if st.lstrip().startswith(filetime_ramsin_varname)][0].lstrip(f"{filetime_ramsin_varname} =").split(".,")[ngrid-1])

        tdict = {"dt": dtlong, "afiletime": afiletime}
#                 print(f'RAMS Z-levels: {i} {zm[i]:.1f} {zt[i]:.1f}')
    return tdict
def multbyfactor(vardict, timedict):
    timestep = timedict["dt"]; afiletime = timedict["afiletime"]
    factordict = {"pertimestep": 1/timestep, #Change DTheta into DTheta/Dt
                "gperkg": 1000, #Convert kg/kg to g/kg
                "mgperkgpersec": 10**6/afiletime, #Convert process budget rates from kg/kg (summed over analysis file interval) to mg/kg/sec
                "ugperkg": 10**9, #Convert kg/kg to ug/kg
                "preciprate": 3600, #Convert mm/s to mm/hr
                "defaultmomentumbudget": 1/(2*timestep)} #Convert from Dw to Dw/Dt
    for vkey in vardict.keys():
        if type(vardict[vkey].unitfactor) == str:
            vardict[vkey].data = vardict[vkey].data*factordict[vardict[vkey].unitfactor]
        elif type(vardict[vkey].unitfactor) != str:
            vardict[vkey].data = vardict[vkey].data*vardict[vkey].unitfactor
    return vardict
    #This multiplies data by a certain factor right before output (previously, and for all derived variable calculatins, they are stored in SI units). For example, this can convert condensate from kg/kg to g/kg, or convert instantaneous latent heating from K to K/s

def roundvars(vardict):
    for vkey in vardict.keys():
        vardict[vkey].data = vardict[vkey].data.round(decimals = vardict[vkey].decnum)
    return vardict
    #Round data to certain number of decimals. This is done right before output, after multbyfactor. For all derived variable calculations, full precision is used.

def get_cartds(rnameflag, vardict, ccoords):
    mvars = xr.Dataset(coords = ccoords)
    varkeys2d = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2d"]
    varkeyspatch = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2dPatch"]
    varkeyssoil = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2dSoil"]
    varkeyssnow = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2dSnow"]
    varkeys3d = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "3d"]
    for vkey in varkeys3d:
        #rnameflag controls whether "RAMS-style" names should be used (rnameflag = 1) or "verbose-style" names should be used (rnameflag = 0) when building the NetCDF
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = ccoords.drop_dims(["patch", "soillev", "snowlev"]).coords, dims = ["z", "y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
            #The reason for this nonsense of going from a Coordinates object to a Dataset then back to Coordinates is that the Coordinates object currently has no "drop_dims" method like a Dataset does :(
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = ccoords.drop_dims(["patch", "soillev", "snowlev"]).coords, dims = ["z", "y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey] #Get rid of data we've already passed into the dataset to avoid hogging memory
    for vkey in varkeys2d:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = ccoords.drop_dims(["z", "patch", "soillev", "snowlev"]).coords, dims = ["y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = ccoords.drop_dims(["z", "patch", "soillev", "snowlev"]).coords, dims = ["y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey]
    for vkey in varkeyspatch:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = ccoords.drop_dims(["z", "soillev", "snowlev"]).coords, dims = ["patch", "y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = ccoords.drop_dims(["z", "soillev", "snowlev"]).coords, dims = ["patch", "y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey]
    for vkey in varkeyssoil:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = ccoords.drop_dims(["z", "snowlev"]).coords, dims = ["patch", "soillev", "y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = ccoords.drop_dims(["z", "snowlev"]).coords, dims = ["patch", "soillev", "y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey]
    for vkey in varkeyssnow:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = ccoords.drop_dims(["z", "soillev"]).coords, dims = ["patch", "snowlev", "y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = ccoords.drop_dims(["z", "soillev"]).coords, dims = ["patch", "snowlev", "y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey]

    return mvars

def get_presds(rnameflag, vardict, pcoords):
    mvars = xr.Dataset(coords = pcoords)
    varkeys2d = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2d"]
    varkeyspatch = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2dPatch"]
    varkeyssoil = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2dSoil"]
    varkeyssnow = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2dSnow"]
    varkeys3d = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "3d"]
    for vkey in varkeys3d:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = pcoords.drop_dims(["patch", "soillev", "snowlev"]).coords, dims = ["pressure_level", "y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
            #The reason for this nonsense of going from a Coordinates object to a Dataset then back to Coordinates is that the Coordinates object currently has no "drop_dims" method like a Dataset does :(
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = pcoords.drop_dims(["patch", "soillev", "snowlev"]).coords, dims = ["pressure_level", "y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey] #Get rid of data we've already passed into the dataset to avoid hogging memory
    for vkey in varkeys2d:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = pcoords.drop_dims(["pressure_level", "patch", "soillev", "snowlev"]).coords, dims = ["y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = pcoords.drop_dims(["pressure_level", "patch", "soillev", "snowlev"]).coords, dims = ["y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey]
    for vkey in varkeyspatch:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = pcoords.drop_dims(["pressure_level", "soillev", "snowlev"]).coords, dims = ["patch", "y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = pcoords.drop_dims(["pressure_level", "soillev", "snowlev"]).coords, dims = ["patch", "y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey]
    for vkey in varkeyssoil:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = pcoords.drop_dims(["pressure_level", "snowlev"]).coords, dims = ["patch", "soillev", "y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = pcoords.drop_dims(["pressure_level", "snowlev"]).coord, dims = ["patch", "soillev", "y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey]
    for vkey in varkeyssnow:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = pcoords.drop_dims(["pressure_level", "soillev"]).coords, dims = ["patch", "snowlev", "y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = pcoords.drop_dims(["pressure_level", "soillev"]).coords, dims = ["patch", "snowlev", "y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey]

    return mvars

def get_sigds(rnameflag, vardict, scoords):
    mvars = xr.Dataset(coords = scoords)
    varkeys2d = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2d"]
    varkeyspatch = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2dPatch"]
    varkeyssoil = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2dSoil"]
    varkeyssnow = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2dSnow"]
    varkeys3d = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "3d"]
    for vkey in varkeys3d:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = scoords.drop_dims(["patch", "soillev", "snowlev"]).coords, dims = ["model_level", "y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
            #The reason for this nonsense of going from a Coordinates object to a Dataset then back to Coordinates is that the Coordinates object currently has no "drop_dims" method like a Dataset does :(
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = scoords.drop_dims(["patch", "soillev", "snowlev"]).coords, dims = ["model_level", "y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey] #Get rid of data we've already passed into the dataset to avoid hogging memory
    for vkey in varkeys2d:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = scoords.drop_dims(["model_level", "patch", "soillev", "snowlev"]).coords, dims = ["y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = scoords.drop_dims(["model_level", "patch", "soillev", "snowlev"]).coords, dims = ["y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey]
    for vkey in varkeyspatch:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = scoords.drop_dims(["model_level", "soillev", "snowlev"]).coords, dims = ["patch", "y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = scoords.drop_dims(["model_level", "soillev", "snowlev"]).coords, dims = ["patch", "y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey]
    for vkey in varkeyssoil:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = scoords.drop_dims(["model_level", "snowlev"]).coords, dims = ["patch", "soillev", "y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = scoords.drop_dims(["model_level", "snowlev"]).coords, dims = ["patch", "soillev", "y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey]
    for vkey in varkeyssnow:
        if rnameflag == 0:
            mvars[vardict[vkey].varname] = xr.DataArray(vardict[vkey].data, coords = scoords.drop_dims(["model_level", "soillev"]).coords, dims = ["patch", "snowlev", "y", "x"])
            mvars[vardict[vkey].varname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].varname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        elif rnameflag == 1:
            mvars[vardict[vkey].ramsname] = xr.DataArray(vardict[vkey].data, coords = scoords.drop_dims(["model_level", "soillev"]).coords, dims = ["patch", "snowlev", "y", "x"])
            mvars[vardict[vkey].ramsname].attrs = {"longname": vardict[vkey].longname, "stdname": vardict[vkey].stdname if vardict[vkey].stdname is not None else "None", "units": vardict[vkey].units, "offline_calculated_variable": str(vardict[vkey].offline)}
            if vardict[vkey].offline == True:
                mvars[vardict[vkey].ramsname].attrs["native_variables_used_to_calculate_derived_quantity"] =  vardict[vkey].invar
        del vardict[vkey]

    return mvars


def cartinterp(gridprops, top2d, atop, rawfile, vardict):
    varkeys3d = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "3d" and vardict[vkey].vector == False] #Separate the vector (u,v,w) and scalar quantities, since the vector quantities have special interpolation requirements
    rawzfull = get_zact(gridprops, top2d)
    ccoords = gen_coords_cart(gridprops, top2d, atop)
    if "WC" in vardict.keys():
        print(f"Interpolating WC --> {[vardict['WC'].varname, vardict['WC'].ramsname][rnameflag]}")
        vardict["WC"].data = get_cartw(rawfile["WC"].values, gridprops, top2d, ccoords["z"].values)
    if "WP" in vardict.keys():
        print(f"Interpolating WP --> {[vardict['WP'].varname, vardict['WP'].ramsname][rnameflag]}")
        vardict["WP"].data = get_cartw(rawfile["WP"].values, gridprops, top2d, ccoords["z"].values)
    if "UC" in vardict.keys():
        print(f"Interpolating UC --> {[vardict['UC'].varname, vardict['UC'].ramsname][rnameflag]}")
        vardict["UC"].data = get_cartu(rawfile["UC"].values, gridprops, rawzfull, ccoords["z"].values)
    if "UP" in vardict.keys():
        print(f"Interpolating UP --> {[vardict['UP'].varname, vardict['UP'].ramsname][rnameflag]}")
        vardict["UP"].data = get_cartu(rawfile["UP"].values, gridprops, rawzfull, ccoords["z"].values)
    if "VC" in vardict.keys():
        print(f"Interpolating VC --> {[vardict['VC'].varname, vardict['VC'].ramsname][rnameflag]}")
        vardict["VC"].data = get_cartv(rawfile["VC"].values, gridprops, rawzfull, ccoords["z"].values)
    if "VP" in vardict.keys():
        print(f"Interpolating VP --> {[vardict['VP'].varname, vardict['VP'].ramsname][rnameflag]}")
        vardict["VP"].data = get_cartv(rawfile["VP"].values, gridprops, rawzfull, ccoords["z"].values)

    vdictlist = [{"data": rawfile[ramsname].values, "ramsname": vardict[ramsname].ramsname, "varname": vardict[ramsname].varname} for ramsname in varkeys3d]
    cartsclpart = partial(get_cartscl_multiprocess, rnameflag, rawzfull, ccoords["z"].values)
    varppool = ProcessPoolExecutor(max_workers = 4) #Parallellize the interpolation of scalar quantities
    try:
        interpscllist = list(varppool.map(cartsclpart, vdictlist))
        varppool.shutdown()
    except:
        print("Oh No!")
        varppool.shutdown()

    print("Done interpolation!")
    for i, entry in enumerate(interpscllist):
        vardict[entry["ramsname"]].data = entry["data"] #Interpolated data is always stored in the post-processing code using "RAMS-style" names. It will still be output using "verbose-style" names if the user chooses that option
        interpscllist[i]["data"] = 1 #This is done to prevent data duplication and hogging memory
    del interpscllist
    del vdictlist
    return vardict

def presinterp(gridprops, top2d, atop, rawfile, vardict):
    varkeys3d = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "3d" and vardict[vkey].vector == False]
    rawpresfull = (rawfile["PI"].values/1004)**(1004/287)
    pcoords = gen_coords_pressure(gridprops, top2d, rawpresfull, userpreslvs)
    if "WC" in vardict.keys():
        print(f"Interpolating WC --> {[vardict['WC'].varname, vardict['WC'].ramsname][rnameflag]}")
        vardict["WC"].data = get_presw(rawfile["WC"].values, gridprops, top2d, rawpresfull, pcoords["pressure_level"].values)
    if "WP" in vardict.keys():
        print(f"Interpolating WP --> {[vardict['WP'].varname, vardict['WP'].ramsname][rnameflag]}")
        vardict["WP"].data = get_presw(rawfile["WP"].values, gridprops, top2d, rawpresfull, pcoords["pressure_level"].values)
    if "UC" in vardict.keys():
        print(f"Interpolating UC --> {[vardict['UC'].varname, vardict['UC'].ramsname][rnameflag]}")
        vardict["UC"].data = get_presu(rawfile["UC"].values, gridprops, rawpresfull, pcoords["pressure_level"].values)
    if "UP" in vardict.keys():
        print(f"Interpolating UP --> {[vardict['UP'].varname, vardict['UP'].ramsname][rnameflag]}")
        vardict["UP"].data = get_presu(rawfile["UP"].values, gridprops, rawpresfull, pcoords["pressure_level"].values)
    if "VC" in vardict.keys():
        print(f"Interpolating VC --> {[vardict['VC'].varname, vardict['VC'].ramsname][rnameflag]}")
        vardict["VC"].data = get_presv(rawfile["VC"].values, gridprops, rawpresfull, pcoords["pressure_level"].values)
    if "VP" in vardict.keys():
        print(f"Interpolating VP --> {[vardict['VP'].varname, vardict['VP'].ramsname][rnameflag]}")
        vardict["VP"].data = get_presv(rawfile["VP"].values, gridprops, rawpresfull, pcoords["pressure_level"].values)

    # print(f"Interpolating 3D Scalars at {ftime.strftime('%Y-%m-%d-%H%M%S')}")
    vdictlist = [{"data": rawfile[ramsname].values, "ramsname": vardict[ramsname].ramsname, "varname": vardict[ramsname].varname} for ramsname in varkeys3d]
    pressclpart = partial(get_presscl_multiprocess, rnameflag, rawpresfull, pcoords["pressure_level"].values)
    varppool = ProcessPoolExecutor(max_workers = 4)
    try:
        interpscllist = list(varppool.map(pressclpart, vdictlist))
        varppool.shutdown()
    except:
        print("Oh No!")
        varppool.shutdown()
    print("Done interpolation!")
    for i, entry in enumerate(interpscllist):
        vardict[entry["ramsname"]].data = entry["data"]
        interpscllist[i]["data"] = 1 #This step is done to prevent data duplication and hogging memory
    del interpscllist
    del vdictlist
    return vardict

def siginterp(gridprops, top2d, atop, rawfile, vardict):
    varkeys3d = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "3d" and vardict[vkey].vector == False]
    scoords = gen_coords_sigma(gridprops, top2d, atop)
    if "WC" in vardict.keys():
        print(f"Interpolating WC --> {[vardict['WC'].varname, vardict['WC'].ramsname][rnameflag]}")
        vardict["WC"].data = get_sigw(rawfile["WC"].values, gridprops, top2d, scoords)
    if "WP" in vardict.keys():
        print(f"Interpolating WP --> {[vardict['WP'].varname, vardict['WP'].ramsname][rnameflag]}")
        vardict["WP"].data = get_sigw(rawfile["WP"].values, gridprops, top2d, scoords)
    if "UC" in vardict.keys():
        print(f"Interpolating UC --> {[vardict['UC'].varname, vardict['UC'].ramsname][rnameflag]}")
        vardict["UC"].data = get_sigu(rawfile["UC"].values, gridprops, scoords)
    if "UP" in vardict.keys():
        print(f"Interpolating UP --> {[vardict['UP'].varname, vardict['UP'].ramsname][rnameflag]}")
        vardict["UP"].data = get_sigu(rawfile["UP"].values, gridprops, scoords)
    if "VC" in vardict.keys():
        print(f"Interpolating VC --> {[vardict['VC'].varname, vardict['VC'].ramsname][rnameflag]}")
        vardict["VC"].data = get_sigv(rawfile["VC"].values, gridprops, scoords)
    if "VP" in vardict.keys():
        print(f"Interpolating VP --> {[vardict['VP'].varname, vardict['VP'].ramsname][rnameflag]}")
        vardict["VP"].data = get_sigv(rawfile["VP"].values, gridprops, scoords)

    # print(f"Subsetting 3D Scalars at {ftime.strftime('%Y-%m-%d-%H%M%S')}")
    vdictlist = [{"data": rawfile[ramsname].values, "ramsname": vardict[ramsname].ramsname, "varname": vardict[ramsname].varname} for ramsname in varkeys3d]
    sigsclpart = partial(get_sigscl_multiprocess, rnameflag, scoords["model_level"].values)
    varppool = ProcessPoolExecutor(max_workers = 4)
    try:
        interpscllist = list(varppool.map(sigsclpart, vdictlist))
        varppool.shutdown()
    except:
        print("Oh No!")
        varppool.shutdown()
    print("Done interpolation!")
    for i, entry in enumerate(interpscllist):
        vardict[entry["ramsname"]].data = entry["data"]
        interpscllist[i]["data"] = 1 #This step is done to prevent data duplication and hogging memory
    del interpscllist
    del vdictlist
    return vardict

def postprocess(prepath, postpath, filetype, gridprops, interptype, atop, userpreslvs, kernname, timedict, window, hydropath, uservars, derivvars, rnameflag, ftime):
    # slcpdict = {"sand": 1465000, "loamy_sand": 1407000, "sandy_loam": 1344000,
    #             "silt_loam": 1273000, "loam": 1214000, "sandy_clay_loam": 1177000,
    #             "silty_clay_loam": 1319000, "clay_loam": 1227000, "sandy_clay": 1177000,
    #             "silty_clay": 1151000, "clay": 1088000, "peat": 874000};
    #This is for a "soil temperature" variable, which I still need to add.
    # afiletime = timedict["afiletime"]; timestep = timedict["dt"]

    st = perf_counter()
    print(f"Post-Processing data at {ftime.strftime('%Y-%m-%d-%H%M%S')}")
    print(f"Reading File {prepath}a-{filetype}-{ftime.strftime('%Y-%m-%d-%H%M%S')}-g{ngrid}.h5")
    rawfile = xr.open_dataset(f"{prepath}a-{filetype}-{ftime.strftime('%Y-%m-%d-%H%M%S')}-g{ngrid}.h5", engine = "h5netcdf", phony_dims = "sort")
    #Get xarray dataset of RAMS output
    top2d = rawfile["TOPT"].values
    datavarlist = list(rawfile.data_vars)
    gridprops["glat"] = rawfile["GLAT"].values; gridprops["glon"] = rawfile["GLON"].values #Assign 2D arrays of latitude and longitude to gridprops
    if uservars == "all":
        uservars = datavarlist

    fullvarlist = fillvarlist()
    vardict = gen_vardict(uservars, datavarlist, fullvarlist); del fullvarlist

    #!!! 
    #Very important! vardict is where all the data for variables (native RAMS variables and derived variables) is stored in the post-processing code. This a dictionary of the form
    #{"rams-style name": outvar}, where outvar is a custom class containing data for the variable, the variable's RAMS-style name, the verbose-style name, the longname, a standard CF-compliant name (if applicable), the variable type (3D, 2D, 2DPatch, 2DSoil, or 2DSnow), the factor by which to multiply the variable before output, the units of the variable, if the varible is a vector, and which variables are used to calculate it (if it's a derived variable). All variables are stored in vardict using their *"RAMS-style" name*, not their *"verbose-style" name*!!
    #!!!

    varkeys2d = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2d"]
    varkeyspatch = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2dPatch"]
    varkeyssoil = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2dSoil"]
    varkeyssnow = [vkey for vkey in vardict.keys() if vardict[vkey].vartype == "2dSnow"]
    interptype = interptype.lower()
    if interptype == "cart":
        ccoords = gen_coords_cart(gridprops, top2d, atop)
        print("Interpolating 3D Scalars to Cartesian Coordinates")
        vardict = cartinterp(gridprops, top2d, atop, rawfile, vardict)
        
    elif interptype == "pressure":
        rawpresfull = (rawfile["PI"].values/1004)**(1004/287)
        pcoords = gen_coords_pressure(gridprops, top2d, rawpresfull, userpreslvs)
        print("Interpolating 3D Scalars to Pressure Coordinates")
        vardict = presinterp(gridprops, top2d, atop, rawfile, vardict)

    elif interptype == "sigma":
        scoords = gen_coords_sigma(gridprops, top2d, atop)
        print("Subsetting 3D Scalars on Sigma-Z Grid")
        vardict = siginterp(gridprops, top2d, atop, rawfile, vardict)

    print("Finished Interpolation!")
    et = perf_counter()
    print(f"Interpolation of variables at time {ftime.strftime('%Y-%m-%d-%H%M%S')} took {et-st:.2f} seconds to calculate")
    for var2d in varkeys2d:
        print(f"Doing variable {var2d}")
        vardict[var2d].data = rawfile[var2d][1:-1,1:-1].values
    for varpatch in varkeyspatch:
        print(f"Doing variable {varpatch}")
        vardict[varpatch].data = rawfile[varpatch][:,1:-1,1:-1].values
    for varsoil in varkeyssoil:
        print(f"Doing Variable {varsoil}")
        vardict[varsoil].data = rawfile[varsoil][:,:,1:-1,1:-1].values
    for varsnow in varkeyssnow:
        print(f"Doing variable {varsnow}")
        vardict[varsnow].data = rawfile[varsnow][:,:,1:-1,1:-1].values

    print("Finished Getting Surface Variables")
    if derivvars is not None:
        print("Finding Derived Quantities")
        if interptype == "cart":
            vardict = get_derivedvars(vardict, derivvars, rawfile, gridprops, window, kernname, rnameflag, hydropath, ccoords = ccoords)
        else:
            vardict = get_derivedvars(vardict, derivvars, rawfile, gridprops, window, kernname, rnameflag, hydropath)

    st = perf_counter()
    print("Beginning Unit Conversion")
    vardict = multbyfactor(vardict, timedict) #Multiply units by the factor specified in the "unitfactor" field of their "outvar" class, to convert from kg/kg to g/kg or K to K/S, etc.
    et = perf_counter()
    print(f"Unit Conversion at {ftime.strftime('%Y-%m-%d-%H%M%S')} took {et-st:.2f} seconds to calculate")
    print("Rounding variables to save disk space")
    vardict = roundvars(vardict) #Round variables
    if interptype == "cart":
        mvars = get_cartds(rnameflag, vardict, ccoords) #Build xarray dataset to save as a NetCDF. Here, variables will be assigned using either "RAMS-style" names or "verbose-style" names, depending on the user's choice. This choice is reflected in "rnameflag". "rnameflag" is 0 if the user chooses "verbose-style" names, and 1 if the user chooses "RAMS-style" names.
    elif interptype == "pressure":
        mvars = get_presds(rnameflag, vardict, pcoords)
    elif interptype == "sigma":
        mvars = get_sigds(rnameflag, vardict, scoords)
    print(f"Saving file at time {ftime.strftime('%Y-%m-%d-%H%M%S')}!")
    comp = {"compression": "gzip", "compression_opts": 6} #Compress the file slightly to save space
    xrenc = {var: comp for var in mvars.data_vars}
    if not os.path.exists(f"{postpath}"): #Make a directory for the post-processed files if one doesn't already exist
        os.mkdir(f"{postpath}")
    mvars.to_netcdf(f"{postpath}mvars-{interptype}-{ftime.strftime('%Y-%m-%d-%H%M%S')}-g{ngrid}.nc", engine = "h5netcdf", encoding = xrenc) #Save files as .nc, not .h5
    mvars.close(); del mvars
    print(f"File at {ftime.strftime('%Y-%m-%d-%H%M%S')} saved!")
    return ftime.strftime('%Y-%m-%d-%H%M%S')
# In[21]:
# if __name__ == "__main__":

starttime = perf_counter()
isuserfile = input("Is there a namelist file you'd like to use for post-processing settings? Yes or No? ")
if isuserfile.strip().upper() not in ["YES", "NO"]:
    raise Exception("Please answer 'yes' or 'no'. ")

if isuserfile.strip().upper() == "YES":
    userfpath = input("Enter the path of the post-processing settings file: ")
    userfpath = userfpath.strip()
    userdict = read_userfile(userfpath)
    folderpath = userdict["folderpath"]
    postpath = userdict["postpath"]
    ramsinpath = userdict["ramsinpath"]
    interptype = userdict["interptype"]
    atop = userdict["atop"]
    userpreslvs = userdict["userpreslvs"]
    instime = userdict["instime"]
    inetime = userdict["inetime"]
    numworkers = int(userdict["numworkers"])
    ngrid = int(userdict["ngrid"])
    uservarin = userdict["uservars"]
    rnameflag = userdict["rnameflag"]
    derivvarin = userdict["derivvars"]
    ywininput = userdict["ywininput"]
    xwininput = userdict["xwininput"]
    kernname = userdict["kernname"]
    filetype = userdict["filetype"]

elif isuserfile.strip().upper() == "NO":
    folderpath = input("Enter the folder name containing simulation output you'd like to analyze: ")
    postpath = input("Enter the folder name where you'd like to save the processed data (the folder will be created if it doesn't currently exist): ")
    ramsinpath = input("Enter the path to the RAMSIN used for this simulation: ")
    interptype = input("Enter the type of vertical coordinates to use: Cartesian coordinates, Pressure coordinates, or Sigma coordinates ")
    instime = input("Enter the analysis start time in YYYY-mm-dd-HHMMSS format: ")
    inetime = input("Enter the analysis end time in YYYY-mm-dd-HHMMSS format: ")
    numworkers = int(input("Number of cores to use for post-processing"))
    ngrid = int(input("Enter the Grid Number to post-process: ")) #For now, can only do one grid at a time; if you want to analyze a fine and coarse grid, will need to run this twice with different ngrid
    uservarin = input("Enter the list of RAMS variables you'd like to process, as comma-separated entries, or put *all* to process all available variables: ")
    rnameflag = input("Use the *RAMS* variable names, or *verbose* variable names in post-processed NetCDF files? ")
    derivvarin = input("Enter the list of derived variables you'd like to output, as comma-separated entries (a full list of available derived variables is avilable in comments at the top of derivedvars.py). Leave blank to not output any derived quantities. Enter *all* to output all derived quantities available for your output file. Enter *nomomentum* to output all variables except momentum budgets, which are quite slow to calculate: ")
if not os.path.exists(folderpath):
    raise Exception("Parent folder of output file directory not found!")

if not os.path.exists(ramsinpath):
    raise Exception("RAMSIN not found!")

if interptype.lower().replace("coordinates", "").rstrip() in ["cartesian", "cart"]:
    interptype = "cart"
elif interptype.lower().replace("coordinates", "").rstrip() in["pressure", "pres"]:
    interptype = "pressure"
elif interptype.lower().replace("coordinates", "").rstrip() in ["sigma", "sig"]:
    interptype = "sigma"
else:
    raise Exception("Please enter a valid interpolation type: Cartesian, Pressure, or Sigma")
if interptype == "cart" or interptype == "sigma":
    if isuserfile.strip().upper() == "NO":
        atop = input("For sigma or cartesian interpolation, please enter the maximum altitude (in meters) you'd like to include in the post-processing. Leave this blank to use the entire grid: ")
    if atop == "":
        atop = None
    else:
        if atop.isnumeric():
            atop = int(atop)
        else:
            raise Exception("atop must be a number or blank!")
    userpreslvs = None
elif interptype == "pressure":
    if isuserfile.strip().upper() == "NO":
        userpreslvs = input("For pressure interpolation, please enter the pressure levels in hPa you'd like the post-processed data on as comma-separated values: ")
    userpreslvs = sorted([int(preslv.strip()) for preslv in userpreslvs.split(",")], reverse = True)
t0 = datetime.strptime(instime, "%Y-%m-%d-%H%M%S")
tf = datetime.strptime(inetime, "%Y-%m-%d-%H%M%S")
if tf<=t0:
    raise Exception("Analysis end time must be later than analysis start time!")
headpath = f"{folderpath}/a-A-{t0.strftime('%Y-%m-%d-%H%M%S')}-head.txt"
prepath = f"{folderpath}/"
postpath = f"{postpath}/"
ftype = filetype.strip("* ").upper()
gridprops = get_gprops(headpath, ngrid)
timedict = get_timedict(ramsinpath, ngrid, ftype)
if uservarin.strip("* ") != "all":
    uservars = [st.strip().upper() for st in uservarin.split(",")]
else:
    uservars = "all"
if rnameflag.strip().upper() == "RAMS":
    rnameflag = 1
elif rnameflag.strip().upper() == "VERBOSE":
    rnameflag = 0
if derivvarin == "":
    derivvars = None
elif derivvarin.strip("* ").lower() == "all":
    dvfile = open("dvardictfile", "rb")
    dvarnamedict = pickle.load(dvfile)
    dvfile.close()
    derivvars = list(dvarnamedict.values())
elif derivvarin.strip("* ").lower() == "nomomentum":
    dvfile = open("dvardictfile", "rb")
    dvarnamedict = pickle.load(dvfile)
    dvfile.close()
    derivvars = list(dvarnamedict.values())
    for momvar in ["RHOP", "BUOY_RHO", "BUOY_THETA", "BUOY_TEMP", "BUOY_VAP", "BUOY_COND", "BUOY_PPRIME", "VPPGF"]:
        derivvars.remove(momvar)

elif derivvarin is not None:
    dvfile = open("dvardictfile", "rb")
    dvarnamedict = pickle.load(dvfile)
    dvfile.close()
    if rnameflag == 0:
        derivvars = [dvarnamedict[st.strip()] for st in derivvarin.split(",")] #This converts from the "verbose-style" derived variable names that will be output to the "RAMS-style" variable names used in the backend
    elif rnameflag == 1:
        derivvars = [st.strip().upper() for st in derivvarin.split(",")]

if isuserfile.strip().upper() == "NO" and any(momvar in derivvars for momvar in ["RHOP", "BUOY_RHO", "BUOY_THETA", "BUOY_TEMP", "BUOY_VAP", "BUOY_COND", "BUOY_PPRIME", "VPPGF"]):
    ywininput = input("How many Y grid points do you want to use for horizontal averaging (for momentum budgets)? ")
    xwininput = input("How many X grid points do you want to use for horizontal averaging (for momentum budgets)? ")
if any(momvar in derivvars for momvar in ["RHOP", "BUOY_RHO", "BUOY_THETA", "BUOY_TEMP", "BUOY_VAP", "BUOY_COND", "BUOY_PPRIME", "VPPGF"]):
    ywindow = int(ywininput.strip())
    xwindow = int(xwininput.strip())
    window = {"ywindow": ywindow, "xwindow": xwindow} #Number of gridpoints to average over in zonal and meridional directions, respectively, when making reference states for momentum budget calculations. 
else:
    window = None

if isuserfile.strip().upper() == "NO":
    kernname = "trikernel" #Change this to "domekernel", "gausskernel", "flatkernel", or "hornkernel" depending on which one you want. Kernel details are described derivcalcs.get_rollvars
hydropath = "hydroparams.csv"
tlist = date_range(instime, inetime, freq = timedelta(seconds = timedict["afiletime"])).to_pydatetime()
postpartial = partial(postprocess, prepath, postpath, ftype, gridprops, interptype, atop, userpreslvs, kernname, timedict, window, hydropath, uservars, derivvars, rnameflag)
with ProcessPoolExecutor(max_workers = numworkers) as ppool: #Parallellize the post-processing routine across time
    flist = print(list(ppool.map(postpartial, tlist)))
print("Running Post-Processing Routine")
# for time in tlist:
#     pst = perf_counter()
#     ft = postpartial(time)
#     pet = perf_counter()
#     print(f"Post-Processing data at time {time.strftime('%Y-%m-%d-%H%M%S')} took {pet-pst:.2f} seconds")

endtime = perf_counter()
print(f"All Post-Processing took {endtime-starttime:.2f} seconds")