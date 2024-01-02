#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import numpy as np
from time import perf_counter
import warnings

def get_gprops(headpath, ngrid):
    # Get grid properties from RAMS header file
    with open(headpath) as f:
        mylist = f.read().splitlines()
        nx = int(mylist[mylist.index('__nnxp')+ngrid+1]) #Number of x grid points on the grid
        ny = int(mylist[mylist.index('__nnyp')+ngrid+1]) #Number of y grid points on the grid
        nz = int(mylist[mylist.index('__nnzp')+ngrid+1]) #Number of z levels on the grid
        dx = float(mylist[mylist.index('__deltaxn')+ngrid+1]) #Grid spacing on the grid
        # dz0 = float(mylist[mylist.index('__deltazn')+ngrid+1]) #Not sure if we need this for anything
        zmix = mylist.index(f'__zmn{str(ngrid).zfill(2)}') 
        zmlvs = np.asarray([float(i) for i in mylist[zmix+2:zmix+nz+2]]) #Momentum (w) z-levels on the grid
        ztix = mylist.index(f'__ztn{str(ngrid).zfill(2)}')
        ztlvs = np.asarray([float(i) for i in mylist[ztix+2:ztix+nz+2]]) #Scalar z-levels on the grid
        dzmlvs = np.diff(zmlvs); dztlvs = np.diff(ztlvs) #Find cell vertical thicknesses
        nsoil = int(mylist[mylist.index('__nzg')+2]) #Number of soil levels
        zsix = mylist.index('__slz') 
        zsoil = np.asarray([float(i) for i in mylist[zsix+2:zsix+nsoil+2]]) #Depth of soil levels
        nsnow = int(mylist[mylist.index('__nzs')+2]) #Number of snow levels
        npatch = int(mylist[mylist.index('__npatch')+2]) #Number of patches in each grid cell

        gridprops = {"nx": nx, "ny": ny, "nz": nz, "dx": dx, "zstarlvs": ztlvs, "dztn": dztlvs, "wzstarlvs": zmlvs, "dzmn": dzmlvs, "nsoil": nsoil, "soillevs": zsoil, "nsnow": nsnow, "npatch": npatch}
    return gridprops

def get_zact(gridprops, top2d): #This returns the altitude of every point on the sigma-z grid
    wzstarlvs = gridprops["wzstarlvs"] #Standard w z-levels
    zstarlvs = gridprops["zstarlvs"] #Standard scalar z-levels
    zmax = wzstarlvs[-1] #Model top altitude
    zstar3d = zstarlvs[:, None, None] #Broadcast standard z levels to 3D
    top3d = np.broadcast_to(top2d, (len(zstarlvs), top2d.shape[0], top2d.shape[1])) #Broadcast topography height to 3D
    zact = zstar3d+top3d*(1-zstar3d/zmax) #Calculate true height of scalar z-levels on sigma grid
    return zact #This returns the FULL array of sigma-z heights, including model bottom, top and sides

def get_newzstar(gridprops, top2d, atop = None):
    #atop is the height we want to clip to, NOT the model top height
    dzms = gridprops["dzmn"]
    if atop is None:
        warnings.warn(f"atop is not set! Using model top as upper boundary!") #If atop is not specified, this will return all heights from the surface to the model top
        constzs = np.arange(top2d.min(), top2d.max()+dzms[0], dzms[0])
        stretchzs = gridprops["wzstarlvs"][1:]+constzs[-1]
        newzlvs = np.concatenate([constzs, stretchzs[stretchzs<gridprops["wzstarlvs"][-1]]])
    elif atop<=top2d.max():
        warnings.warn(f"Maxmimum terrain height is {top2d.max():.2f}m, but atop is {atop}m! Data over high terrain will be lost!") #If atop is below the height of maximum topography, some information over high terrain will be lost when using cartesian interpolation
        constzs = np.arange(top2d.min(), atop+1, dzms[0])
        return constzs
    elif atop<=gridprops["wzstarlvs"][-1]:
        constzs = np.arange(top2d.min(), top2d.max()+dzms[0], dzms[0])
        stretchzs = gridprops["wzstarlvs"][1:]+constzs[-1]
        newzlvs = np.concatenate([constzs, stretchzs[stretchzs<atop]])
    else:
        warnings.warn(f"atop is {atop}m, but model top is {gridprops['wzstarlvs'][-1]}m! Truncating analysis to model top!") #If the user attempts to get interpolated output which extends higher than the model top, truncate to the model top
        constzs = np.arange(top2d.min(), top2d.max()+dzms[0], dzms[0])
        stretchzs = gridprops["wzstarlvs"][1:]+constzs[-1]
        newzlvs = np.concatenate([constzs, stretchzs[stretchzs<gridprops["wzstarlvs"][-1]]])
    return newzlvs

def get_zlvs_sigma(gridprops, top2d, atop = None):
    mintopocoords = divmod(np.argmin(top2d), top2d.shape[1]) #Unfortunately, np.argmin returns the index of the minimum in the flattened array, so divmod is necessary to get the 2D index
    zact = get_zact(gridprops, top2d) #3D array of true heights on the sigma-z grid
    mintopozact = zact[:,mintopocoords[0], mintopocoords[1]] #Sigma-z levels over lowest-topography area. If lowest topography is sea level, this will be equivalent to zstarlvs
    wzstarlvs = gridprops["wzstarlvs"]
    if atop is None:
        warnings.warn(f"atop is not defined! Using model top as upper boundary!")
        zactsub = zact[1:,:,:]
        zlvs = np.arange(1, gridprops["nz"])
    elif atop<=top2d.min():
        warnings.warn(f"Minimum terrain height is {top2d.max():.2f}m, but atop is {atop}m! Only surface level will be output!") #If the topography everywhere in the domain is higher than the user-set atop (for example if all terrain is above sea level and the user selects an atop of 0), return only the lowest atmospheric model level
        zactsub = zact[1,:,:] #Zeroth z-level on scalar grid is below terrain; First z-level is dz/2 above terrain
        zlvs = 1
    elif atop<=top2d.max():
        warnings.warn(f"Maximum terrain height is {top2d.max():.2f}m, but atop is {atop}m!")
        zlvs = np.arange(1, np.argmax(mintopozact[mintopozact<=atop])) #Find the model level where the sigma-z height over the minimum topography is equal to atop
        zactsub = zact[zlvs,:,:] #Because the sigma-z grid is terrain-following, the height of the top sigma-z level over the lowest topography will be less than or equal to atop, while the height of the top sigma-z level over the highest topography will be greater than atop
    elif atop<=wzstarlvs[-1]:
        zlvs = np.arange(1, np.argmax(mintopozact[mintopozact<=atop]))
        zactsub = zact[zlvs,:,:]
    else:
        warnings.warn(f"atop is {atop}m, but model top height is {wzstarlvs[-1]}m! Truncating analysis to model top!")
        zactsub = zact[1:,:,:]
        zlvs = np.arange(1, gridprops["nz"])
    return {"zlevels": zlvs, "zactual": zactsub}

def gen_coords_sigma(gridprops, top2d, atop):
    glat = gridprops["glat"]; glon = gridprops["glon"] #2D latitude and longitude arrays
    nx = gridprops["nx"]; ny = gridprops["ny"]; dx = gridprops["dx"] #x and y distance (in meters) from grid center
    zdict = get_zlvs_sigma(gridprops, top2d, atop) #3D array of true altitude along sigma-z levels
    sclxcoords = np.arange(-(np.ceil(nx/2)-1.5)*dx, (np.floor(nx/2)-0.5)*dx, dx) #x-coordinates of scalar variables (in the center of model grid cells)
    sclycoords = np.arange(-(np.ceil(ny/2)-1.5)*dx, (np.floor(ny/2)-0.5)*dx, dx) #y-coordinates of scalar variables (in the center of model grid cells)
    mergedscl = xr.Dataset(coords = {"model_level": zdict["zlevels"], "y": sclycoords, "x": sclxcoords})
    mergedscl["model_level"].attrs = {"longname": "Model Sigma Level", "units": 1}; mergedscl["y"].attrs = {"longname": "Northward Distance from Model Center", "units": "m"}; mergedscl["x"].attrs = {"longname": "Eastward Distance from Model Center", "units": "m"}
    mergedscl = mergedscl.assign_coords({"lat1d": ("y", glat[1:-1,int(nx/2)]), "lon1d": ("x", glon[int(ny/2),1:-1]), "lat2d": (["y", "x"], glat[1:-1, 1:-1]), "lon2d": (["y", "x"], glon[1:-1, 1:-1]), "z": (["model_level", "y", "x"], zdict["zactual"][:,1:-1,1:-1]),
    "soillev": gridprops["soillevs"], "patch": np.arange(0, gridprops["npatch"]), "snowlev": np.arange(0, gridprops["nsnow"])})
    #lat1d and lon1d are approximate latitude and longitude coordinates for if the user wants to assign these as dimension coordinates in an xarray dataset. While the latitude array is not perfectly constant with longitude, and the longitude array is not perfectly constant with latitude, these 1d arrays will be relatively accurate, especially for small domains or those near the equator. lon2d and lat2d are the full 2D arrays of longitude and latitude, if the user prefers to use those instead.

    mergedscl["z"].attrs = {"longname": "Altitude Above Mean Sea Level on Sigma-Z Grid", "units": "m"}
    mergedscl["lat1d"].attrs = {"longname": "1D Latitude (Based on Middle Longitude of the Model Grid)", "units": "degree_north"}; mergedscl["lon1d"].attrs = {"longname": "1D Longitude (Based on Middle Latitudes of the Model Grid)", "units": "degree_east"}; mergedscl["lat2d"].attrs = {"longname": "2D Latitude Field from Full Model Grid", "units": "degree"}; mergedscl["lon2d"].attrs = {"longname": "2D Longitude Field from the Full Model Grid", "units": "degree"}
    mergedscl["soillev"].attrs = {"longname": "Depth of Soil Below Surface", "units": "m"}; mergedscl["patch"].attrs = {"longname": "Patch Number (0 = water, 1+ = land)", "units": 1}; mergedscl["snowlev"].attrs = {"longname": "Snow Level Number", "units": 1}
    #Unlike the raw RAMS output, all post-processed data will be associated with a explicitly defined dimensions and coordinates (no need to worry about what "phony_dim_1" means anymore!). These coordinates and dimensions have units attached, which makes xarray operations such as differentiation and integration easy.
    return mergedscl
    #The vertical dimension of this dataset is model level. There is a 3D coordinate called "z" which gives the true 3D array of terrain-following heights, but because xarray only supports indexing and arithmetic with dimensions and not coordinates, this 3D array of altitude cannot be used for builtin xarray operations such as .sel, .interp, .integrate, or .differentiate
    
def gen_coords_pressure(gridprops, top2d, presvals, userpreslvs):
    glat = gridprops["glat"]; glon = gridprops["glon"]
    nx = gridprops["nx"]; ny = gridprops["ny"]; dx = gridprops["dx"]
    zact = get_zact(gridprops, top2d)[1:-1, 1:-1, 1:-1] #Trim model bottom, top, and lateral boundaries for terrain-following height array
    zstarpreslev = np.zeros((len(userpreslvs), zact.shape[1], zact.shape[2]))
    presvals = presvals[1:-1,1:-1,1:-1] #Trim model bottom, top, and lateral boundaries for raw pressure data.
    for j in range(zact.shape[1]):
        for i in range(zact.shape[2]):
            zstarpreslev[:,j,i] = np.interp(userpreslvs[::-1], presvals[::-1,j,i], zact[::-1,j,i], right = np.nan)[::-1] #Get 3D array of altitude along constant pressure levels
    sclxcoords = np.arange(-(np.ceil(nx/2)-1.5)*dx, (np.floor(nx/2)-0.5)*dx, dx)
    sclycoords = np.arange(-(np.ceil(ny/2)-1.5)*dx, (np.floor(ny/2)-0.5)*dx, dx)
    mergedscl = xr.Dataset(coords = {"pressure_level": userpreslvs, "y": sclycoords, "x": sclxcoords})
    mergedscl["pressure_level"].attrs = {"longname": "Pressure Level", "units": "hPa"}; mergedscl["y"].attrs = {"longname": "Northward Distance from Model Center", "units": "m"}; mergedscl["x"].attrs = {"longname": "Eastward Distance from Model Center", "units": "m"}
    mergedscl = mergedscl.assign_coords({"lat1d": ("y", glat[1:-1,int(nx/2)]), "lon1d": ("x", glon[int(ny/2),1:-1]), "lat2d": (["y", "x"], glat[1:-1, 1:-1]), "lon2d": (["y", "x"], glon[1:-1, 1:-1]), "z": (["pressure_level", "y", "x"], zstarpreslev),
    "soillev": gridprops["soillevs"], "patch": np.arange(0, gridprops["npatch"]), "snowlev": np.arange(0, gridprops["nsnow"])})
    mergedscl["soillev"].attrs = {"longname": "Depth of Soil Below Surface", "units": "m"}; mergedscl["patch"].attrs = {"longname": "Patch Number (0 = water, 1+ = land)", "units": 1}; mergedscl["snowlev"].attrs = {"longname": "Snow Level Number", "units": 1}
    mergedscl["lat1d"].attrs = {"longname": "1D Latitude (Based on Middle Longitude of the Model Grid)", "units": "degree_north"}; mergedscl["lon1d"].attrs = {"longname": "1D Longitude (Based on Middle Latitudes of the Model Grid)", "units": "degree_east"}; mergedscl["lat2d"].attrs = {"longname": "2D Latitude Field from Full Model Grid", "units": "degree"}; mergedscl["lon2d"].attrs = {"longname": "2D Longitude Field from the Full Model Grid", "units": "degree"}
    return mergedscl
    #The vertical dimension of this dataset is pressure level, in hPa. This means that it is easy to select data along a certain pressure level by doing ".sel(pressure_level = 500)" or something like that. As in the sigma-z dataset above, there is a 3D "z" coordinate which gives height above sea level along constant pressure levels.

def gen_coords_cart(gridprops, top2d, atop):
    glat = gridprops["glat"]; glon = gridprops["glon"];
    nx = gridprops["nx"]; ny = gridprops["ny"]; dx = gridprops["dx"]
    newzlvs = get_newzstar(gridprops, top2d, atop)
    sclxcoords = np.arange(-(np.ceil(nx/2)-1.5)*dx, (np.floor(nx/2)-0.5)*dx, dx)
    sclycoords = np.arange(-(np.ceil(ny/2)-1.5)*dx, (np.floor(ny/2)-0.5)*dx, dx)
    mergedscl = xr.Dataset(coords = {"z": newzlvs[1:], "y": sclycoords, "x": sclxcoords})
    mergedscl["z"].attrs = {"longname": "Altutude Above Mean Sea Level", "unitss": "m"}; mergedscl["y"].attrs = {"longname": "Northward Distance from Model Center", "units": "m"}; mergedscl["x"].attrs = {"longname": "Eastward Distance from Model Center", "units": "m"}
    mergedscl = mergedscl.assign_coords({"lat1d": ("y", glat[1:-1,int(nx/2)]), "lon1d": ("x", glon[int(ny/2),1:-1]), "lat2d": (["y", "x"], glat[1:-1, 1:-1]), "lon2d": (["y", "x"], glon[1:-1, 1:-1]),
    "soillev": gridprops["soillevs"], "patch": np.arange(0, gridprops["npatch"]), "snowlev": np.arange(0, gridprops["nsnow"])})
    mergedscl["lat1d"].attrs = {"longname": "1D Latitude (Based on Middle Longitude of the Model Grid)", "units": "degree_north"}; mergedscl["lon1d"].attrs = {"longname": "1D Longitude (Based on Middle Latitudes of the Model Grid)", "units": "degree_east"}; mergedscl["lat2d"].attrs = {"longname": "2D Latitude Field from Full Model Grid", "units": "degree"}; mergedscl["lon2d"].attrs = {"longname": "2D Longitude Field from the Full Model Grid", "units": "degree"}
    mergedscl["soillev"].attrs = {"longname": "Depth of Soil Below Surface", "units": "m"}; mergedscl["patch"].attrs = {"longname": "Patch Number (0 = water, 1+ = land)", "units": 1}; mergedscl["snowlev"].attrs = {"longname": "Snow Level Number", "units": 1}
    return mergedscl
    #The vertical dimension of this dataset is height above sea level, in m. This means that it is easy to select data long a constant height by doing ".sel(z = 1000)" or something like that. This also makes it easy to do vertical integration, differentiation, and interpolation in cartesian space, which cannot be done when using pressure or sigma-z output.


def get_cartw(wvalsraw, gridprops, top2d, newzlvs):
    zstarlvs = gridprops["zstarlvs"] #Sigma-z scalar levels (first level is under the surface)
    wzstarlvs = gridprops["wzstarlvs"] #Sigma-z momentum (w) levels (first level is at 0)
    zmax = wzstarlvs[-1] #Model top height is determined by the highest momentum sigma-z level
    top3d = np.broadcast_to(top2d, (len(zstarlvs), top2d.shape[0], top2d.shape[1]))
    wzstar3d = np.broadcast_to(wzstarlvs[:, None, None], (len(wzstarlvs), top2d.shape[0], top2d.shape[1]))
    #These are just both placeholder arrays for calculating true 3d heights on the momentum z-levels
    wzstar3d = wzstar3d[0:-1, 1:-1, 1:-1] #Truncate to exclude the top levels and the horizontal boundaries
    wzact = wzstar3d+top3d[0:-1, 1:-1, 1:-1]*(1-wzstar3d/zmax) #Slight reformulation compared to (but mathematically equivalent to) the way this is done in REVU
    wvals = wvalsraw[0:-1, 1:-1, 1:-1]
    winterp = np.zeros((len(newzlvs), wvals.shape[1], wvals.shape[2]), dtype = np.float32)
    #For every i and j horizontal point, interpolate the true altitudes of the SCALAR levels to a 1-d set of standard altitudes
    for i in range(0, winterp.shape[2]):
        for j in range(0, winterp.shape[1]):
            winterp[:,j,i] = np.interp(newzlvs, wzact[:,j,i], wvals[:,j,i], left = np.nan) #Notice our interpolation will have more levels than our original data
    #w is defined by z,y,x, but has additional coordinates for lat and lon, in both 1d and 2d.                 
    return winterp
    #This returns w along cartesian z levels

def get_presw(wvalsraw, gridprops, top2d, presvals, userpreslvs):
    zstarlvs = gridprops["zstarlvs"] #Sigma-z scalar levels (first level is under the surface)
    wzstarlvs = gridprops["wzstarlvs"] #Sigma-z momentum (w) levels (first level is at 0)
    zmax = wzstarlvs[-1] #Model top height is determined by the highest momentum sigma-z level
    top3d = np.broadcast_to(top2d, (len(zstarlvs), top2d.shape[0], top2d.shape[1]))
    wzstar3d = np.broadcast_to(wzstarlvs[:, None, None], (len(wzstarlvs), top2d.shape[0], top2d.shape[1]))
    pzact = get_zact(gridprops, top2d)[1:-1, 1:-1, 1:-1] #trim model bottom, top, and lateral boundaries from sigma-z array
    #These are just both placeholder arrays for calculating true 3d heights on the momentum z-levels
    wzstar3d = wzstar3d[0:-1, 1:-1, 1:-1] #Truncate to exclude the top levels and the horizontal boundaries
    wzact = wzstar3d+top3d[0:-1, 1:-1, 1:-1]*(1-wzstar3d/zmax) #Slight reformulation compared to (but mathematically equivalent to) the way this is done in REVU
    wvals = wvalsraw[0:-1, 1:-1, 1:-1]
    presvals = presvals[1:-1,1:-1,1:-1] #Get rid of model top and bottom and lateral boundaries for pressure field
    wpresinterp = np.zeros((wzact.shape[0], wzact.shape[1], wzact.shape[2]), dtype = np.float32)
    winterp = np.zeros((len(userpreslvs), presvals.shape[1], presvals.shape[2]), dtype = np.float32)
    #For every i and j horizontal point, interpolate the true altitudes of the SCALAR levels to a 1-d set of standard altitudes
    for i in range(0, winterp.shape[2]):
        for j in range(0, winterp.shape[1]):
            wpresinterp[:,j,i] = np.interp(wzact[:,j,i], pzact[:,j,i], presvals[:,j,i]) #First, interpolate pressure values from the scalar sigma-z grid to the w sigma-z grid
            winterp[:,j,i] = np.interp(userpreslvs[::-1], wpresinterp[::-1,j,i], wvals[::-1,j,i], right = np.nan)[::-1] #Then, interpolate from these w-grid pressure values to standard pressure levels. We have to do the weird ::-1 thing because np.interp needs index points to increase, and pressure decreases with height. So we need to flip the interpolation upside-down, then flip it back right-side up (the ::-1 on the very end)
    return winterp
    #This returns w along pressure levels

def get_sigw(wvalsraw, gridprops, top2d, scoords):
    zstarlvs = gridprops["zstarlvs"] #Sigma-z scalar levels (first level is under the surface)
    wzstarlvs = gridprops["wzstarlvs"] #Sigma-z momentum (w) levels (first level is at 0)
    zmax = wzstarlvs[-1] #Model top height is determined by the highest momentum sigma-z level
    top3d = np.broadcast_to(top2d, (len(zstarlvs), top2d.shape[0], top2d.shape[1]))
    wzstar3d = np.broadcast_to(wzstarlvs[:, None, None], (len(wzstarlvs), top2d.shape[0], top2d.shape[1]))
    zact = scoords["z"].values #Pull sigma-z heights from scoords, where they have already been computed
    #These are just both placeholder arrays for calculating true 3d heights on the momentum z-levels
    wzstar3d = wzstar3d[0:-1, 1:-1, 1:-1] #Truncate to exclude the top levels and the horizontal boundaries
    wzact = wzstar3d+top3d[0:-1, 1:-1, 1:-1]*(1-wzstar3d/zmax) #Slight reformulation compared to (but mathematically equivalent to) the way this is done in REVU
    wvals = wvalsraw[0:-1, 1:-1, 1:-1]
    winterp = np.zeros((len(scoords["model_level"].values), wvals.shape[1], wvals.shape[2]), dtype = np.float32)
    #For every i and j horizontal point, interpolate the true altitudes of the SCALAR levels to a 1-d set of standard altitudes
    for i in range(0, winterp.shape[2]):
        for j in range(0, winterp.shape[1]):
            winterp[:,j,i] = np.interp(zact[:,j,i], wzact[:,j,i], wvals[:,j,i], left = np.nan) #Notice our interpolation will have more levels than our original data
    # winterpda = xr.DataArray(winterp, coords = scoords.drop_dims(["soillev", "patch", "snowlev"]).coords, name = "w")
    #w is defined by z,y,x, but has additional coordinates for lat and lon, in both 1d and 2d.                 
    return winterp
    #This returns w along model sigma-z levels. This is not exactly the same as the raw RAMS output, because here the w values have been interpolated to the center of the vertical grid cells, rather than being at the top and bottoms of the cells as in the raw output.

def get_cartu(uvalsraw, gridprops, zact, newzlvs):
    nx = gridprops["nx"]; ny = gridprops["ny"]; dx = gridprops["dx"]
    uvals = uvalsraw[1:-1,1:-1,0:-1]
    uinterp_offset = np.zeros((len(newzlvs), uvals.shape[1], uvals.shape[2]), dtype = np.float32) #This will hold u values on cartesian z levels, but still on a staggered horizontal grid
    uinterp_true = np.zeros((len(newzlvs), uvals.shape[1], uvals.shape[2]-1), dtype = np.float32) #This will hold the u values in the centers of the grid cells on the cartesian grid
    sclxcoords = np.arange(-(np.ceil(nx/2)-1.5)*dx, (np.floor(nx/2)-0.5)*dx, dx)
    uzact = zact[1:-1,1:-1,0:-1]
    for i in range(0, uinterp_offset.shape[2]):
        for j in range(0, uinterp_offset.shape[1]):
            uinterp_offset[:,j,i] = np.interp(newzlvs, uzact[:,j,i], uvals[:,j,i], left = np.nan) #Notice our interpolation will have more levels than our original data
    
    #Unlike w, we now need to interpolate in the horizontal because u is on a horizontally-staggered grid. This will interpolate u to the scalar grid points
    for k in range(0, uinterp_offset.shape[0]):
        for j in range(0, uinterp_offset.shape[1]):
            uinterp_true[k,j,:] = np.interp(sclxcoords, np.arange(-(np.ceil(nx/2)-1)*dx, np.floor(nx/2)*dx, dx), uinterp_offset[k,j,:])
    return uinterp_true
    #Return u along z levels, interpolated to the center of horizontal grid cells

def get_presu(uvalsraw, gridprops, presvals, userpreslvs):
    nx = gridprops["nx"]; ny = gridprops["ny"]; dx = gridprops["dx"]
    uvals = uvalsraw[1:-1,1:-1,0:-1]
    presvals = presvals[1:-1, 1:-1, 0:-1]
    uinterp_offset = np.zeros((len(userpreslvs), uvals.shape[1], uvals.shape[2]), dtype = np.float32) #This will hold u values on pressure levels, but still on a staggered horizontal grid
    uinterp_true = np.zeros((len(userpreslvs), uvals.shape[1], uvals.shape[2]-1), dtype = np.float32) #This will hold the u values in the center of grid cells on pressure levels
    for j in range(uvals.shape[1]):
        for i in range(uvals.shape[2]):
            uinterp_offset[:,j,i] = np.interp(userpreslvs[::-1], presvals[::-1,j,i], uvals[::-1,j,i], right = np.nan)[::-1]
   
    #Unlike w, we now need to interpolate in the horizontal because u is on a horizontally-staggered grid. This will interpolate u to the scalar points
    sclxcoords = np.arange(-(np.ceil(nx/2)-1.5)*dx, (np.floor(nx/2)-0.5)*dx, dx)
    for k in range(0, uinterp_offset.shape[0]):
        for j in range(0, uinterp_offset.shape[1]):
            uinterp_true[k,j,:] = np.interp(sclxcoords, np.arange(-(np.ceil(nx/2)-1)*dx, np.floor(nx/2)*dx, dx), uinterp_offset[k,j,:])
    return uinterp_true
    #Return u along constant pressure levels, interpolated to the center of horizontal grid cells

def get_sigu(uvalsraw, gridprops, scoords):
    nx = gridprops["nx"]; ny = gridprops["ny"]; dx = gridprops["dx"]
    uinterp_offset = uvalsraw[1:scoords["model_level"].values[-1]+1,1:-1,0:-1] #U data on the C grid
    uinterp_true = np.zeros((uinterp_offset.shape[0], uinterp_offset.shape[1], uinterp_offset.shape[2]-1), dtype = np.float32) #What will hold U data interpolated to grid point centers
    
    #Unlike w, we now need to interpolate in the horizontal because u is on a horizontally-staggered grid. This will interpolate u to the scalar points
    sclxcoords = np.arange(-(np.ceil(nx/2)-1.5)*dx, (np.floor(nx/2)-0.5)*dx, dx)
    for k in range(0, uinterp_offset.shape[0]):
        for j in range(0, uinterp_offset.shape[1]):
            uinterp_true[k,j,:] = np.interp(sclxcoords, np.arange(-(np.ceil(nx/2)-1)*dx, np.floor(nx/2)*dx, dx), uinterp_offset[k,j,:])
    return uinterp_true
    #Return u along sigma-z levels, interpolated to the center of horizontal grid cells, rather than on the zonal boundaries as in the raw RAMS output

def get_cartv(vvalsraw, gridprops, zact, newzlvs):
    nx = gridprops["nx"]; ny = gridprops["ny"]; dx = gridprops["dx"]
    vvals = vvalsraw[1:-1,0:-1,1:-1]
    vinterp_offset = np.zeros((len(newzlvs), vvals.shape[1], vvals.shape[2]), dtype = np.float32)
    vinterp_true = np.zeros((len(newzlvs), vvals.shape[1]-1, vvals.shape[2]), dtype = np.float32)
    sclycoords = np.arange(-(np.ceil(ny/2)-1.5)*dx, (np.floor(ny/2)-0.5)*dx, dx)
    vzact = zact[1:-1,0:-1,1:-1]
    for i in range(0, vinterp_offset.shape[2]):
        for j in range(0, vinterp_offset.shape[1]):
            vinterp_offset[:,j,i] = np.interp(newzlvs, vzact[:,j,i], vvals[:,j,i], left = np.nan)
    #Unlike w, we need to interpolate in the horizontal because v is on a horizontally-staggered grid. This will interpolate u to the scalar points
    for k in range(0, vinterp_offset.shape[0]):
        for i in range(0, vinterp_offset.shape[2]):
            vinterp_true[k,:,i] = np.interp(sclycoords, np.arange(-(np.ceil(ny/2)-1)*dx, np.floor(ny/2)*dx, dx), vinterp_offset[k,:,i])
    return vinterp_true
    #Same as u above, but for v

def get_presv(vvalsraw, gridprops, presvals, userpreslvs):
    nx = gridprops["nx"]; ny = gridprops["ny"]; dx = gridprops["dx"]
    vvals = vvalsraw[1:-1,0:-1,1:-1]
    presvals = presvals[1:-1, 0:-1, 1:-1]
    vinterp_offset = np.zeros((len(userpreslvs), vvals.shape[1], vvals.shape[2]), dtype = np.float32)
    vinterp_true = np.zeros((len(userpreslvs), vvals.shape[1]-1, vvals.shape[2]), dtype = np.float32)
    for j in range(vvals.shape[1]):
        for i in range(vvals.shape[2]):
            vinterp_offset[:,j,i] = np.interp(userpreslvs[::-1], presvals[::-1,j,i], vvals[::-1,j,i], right = np.nan)[::-1]
    #Unlike w, we need to interpolate in the horizontal because v is on a horizontally-staggered grid. This will interpolate u to the scalar points
    sclycoords = np.arange(-(np.ceil(ny/2)-1.5)*dx, (np.floor(ny/2)-0.5)*dx, dx)
    for k in range(0, vinterp_offset.shape[0]):
        for i in range(0, vinterp_offset.shape[2]):
            vinterp_true[k,:,i] = np.interp(sclycoords, np.arange(-(np.ceil(ny/2)-1)*dx, np.floor(ny/2)*dx, dx), vinterp_offset[k,:,i])
    return vinterp_true
    #Same as u above, but for v

def get_sigv(vvalsraw, gridprops, scoords):
    nx = gridprops["nx"]; ny = gridprops["ny"]; dx = gridprops["dx"]
    vinterp_offset = vvalsraw[1:scoords["model_level"].values[-1]+1,0:-1,1:-1]
    vinterp_true = np.zeros((vinterp_offset.shape[0], vinterp_offset.shape[1]-1, vinterp_offset.shape[2]), dtype = np.float32)
    #Unlike w, we need to interpolate in the horizontal because v is on a horizontally-staggered grid. This will interpolate u to the scalar points
    sclycoords = np.arange(-(np.ceil(ny/2)-1.5)*dx, (np.floor(ny/2)-0.5)*dx, dx)
    for k in range(0, vinterp_offset.shape[0]):
        for i in range(0, vinterp_offset.shape[2]):
            vinterp_true[k,:,i] = np.interp(sclycoords, np.arange(-(np.ceil(ny/2)-1)*dx, np.floor(ny/2)*dx, dx), vinterp_offset[k,:,i])
    return vinterp_true
    #Same as u above, but for v

def get_cartscl_multiprocess(rnameflag, zact, newzlvs, vdict):
    #Simplest of all, since the variables are already on the scalar grid points, not offset like the vector quantities
    sclzact = zact[1:-1,1:-1,1:-1]
    if rnameflag == 0:
        print(f"Interpolating variable {vdict['ramsname']} --> {vdict['varname']}")
    elif rnameflag == 1:
        print(f"Interpolating variable {vdict['ramsname']}")
    sclvals = vdict["data"][1:-1, 1:-1, 1:-1] #Chop off the lateral boundaries
    sclinterp = np.zeros((len(newzlvs), sclvals.shape[1], sclvals.shape[2]), dtype = np.float32)
    for i in range(sclvals.shape[2]):
        for j in range(sclvals.shape[1]):
            sclinterp[:,j,i] = np.interp(newzlvs, sclzact[:,j,i], sclvals[:,j,i], left = np.nan)
    return {"ramsname": vdict["ramsname"], "data": sclinterp} #We just want the values, not a full dataarray. Coordinates will be assigned right at the end
    #Interpolate 3D scalar quantities to constant z levels

def get_presscl_multiprocess(rnameflag, presvals, userpreslvs, vdict):
    if rnameflag == 0:
        print(f"Interpolating variable {vdict['ramsname']} --> {vdict['varname']}")
    elif rnameflag == 1:
        print(f"Interpolating variable {vdict['ramsname']}")
    sclvals = vdict["data"][1:-1,1:-1,1:-1]
    presvals = presvals[1:-1,1:-1,1:-1] #Get rid of model bottom and top and lateral boundaries
    sclinterp = np.zeros((len(userpreslvs), sclvals.shape[1], sclvals.shape[2]))
    for j in range(sclvals.shape[1]):
        for i in range(sclvals.shape[2]):
            sclinterp[:,j,i] = np.interp(userpreslvs[::-1], presvals[::-1,j,i], sclvals[::-1,j,i], right = np.nan)[::-1]
    return {"ramsname": vdict["ramsname"], "data": sclinterp} #We just want the values, not a full dataarray. Coordinates will be assigned right at the end
    #Interpolate 3D scalar quantities to pressure levels

def get_sigscl_multiprocess(rnameflag, levels, vdict):
    if rnameflag == 0:
        print(f"Subsetting variable {vdict['ramsname']} --> {vdict['varname']}")
    elif rnameflag == 1:
        print(f"Subsetting variable {vdict['ramsname']}")
    sclvals = vdict["data"][1:levels[-1]+1,1:-1,1:-1]
    return {"ramsname": vdict["ramsname"], "data": sclvals}
    #Return raw RAMS output on the sigma-z grid, truncating the model bottom, lateral boundaries, and subsetting to the highest sigma-z level requested by the user
    

if __name__ == "__main__":
    # All this is just testing stuff - you can delete it if it's not necessary
    headpath = "/moonbow/ascheb/test.les/2010/hires_control/rams_output/a-A-2010-01-01-120000-head.txt"
    afile = xr.open_dataset("/moonbow/ascheb/test.les/2010/hires_control/rams_output/a-A-2010-01-02-230000-g1.h5", engine = "h5netcdf", phony_dims = "sort")
    print(afile["WC"])
    atop = 10000; ngrid = 1
    gprops = get_gprops(headpath, ngrid)
    gprops["glat"] = afile["GLAT"].values; gprops["glon"] = afile["GLON"].values
    ramstopo = afile["TOPT"]
    winterp = get_cartw(afile["WC"].values, ramstopo.values, gprops, atop)
    print(winterp)
    vinterp = get_cartv(afile["VC"].values, ramstopo.values, gprops, atop)
    print(vinterp)
    thetainterp = get_cartscl(afile, ["THETA", "RV", "PI"], ramstopo.values, gprops, atop)
    print(thetainterp)
    thetainterp.to_netcdf("testtheta.nc", engine = "h5netcdf")
    uinterp = get_cartu(afile["UC"].values, ramstopo.values, gprops, atop)
    print(uinterp)
    mvars_postprocess = xr.merge([winterp, vinterp, uinterp, thetainterp])
    mvars_postprocess.to_netcdf("testmvars.nc", engine = "h5netcdf")