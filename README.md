
=======
# rams_postprocess
Post-Processing code for RAMS model output

Can interpolate variables to cartesian or pressure levels, or keep the data on a sigma-z grid
Interpolates vector quantities (like u, v, w) from grid cell edges of the Arikawa C-grid to the grid cell centers
Returns many derived variables from model output (such as pressure, temperature, dewpoint, particle number concentrations, particle diameters, buoyancy, and more!)
Parallelized using the builtin concurrent.futures package in Python

Required dependencies: python3.9 or later, numpy, pandas, xarray, h5netcdf, astropy

Currently, there's very little documentation. I will work on that going forward. To run the postprocessing code, do the following:
1. Make/modify a namelist file (There's an example file called "postprocessfile" in the directory. Don't change the names of the prompts, but do change the answers on the line below to what you want (where your RAMS output files are, where you want your files to go, what type of interpolation, what variables, etc).
2. Run the postprocessing code with "python3 postprocess.py" in the terminal.
3. When the text prompt asks "Is there a file you'd like to use for post-processing settings? Yes or No?", answer "yes"
4. When prompted for a file path, enter the path to the namelist file you just made.

That's it! Then, the program should run automatically until it has finished post-processing the files for all the times you requested.

Because of the parallelization, I would highly recommend running this code on a server with several cores and a decent amount of memory (for the van den Heever group, that's Snowfall1/2/3 or Solvarg). In the namelist, choose a number of cores you feel is appropriate (I've found that the number of workers you should use to avoid using all up memory is *roughly* (0.6*Total Server Memory in bytes)/(nx*ny*nz*nv*8.5), where nx, ny, nz are the number of x,y, and z points in the post-processed coordinates, nv is the number of 3D variables (which dominate memory use). 

Procedure for adding a new RAMS output variable (*not* a derived variable):
1. Add the new variable in the "get_fullvarlist" function in varinit.py. A guide to the different attributes of the "outvar" class are detailed at the top of varinit.py.

Procedure for adding a new derived variable (*not* a native RAMS output variable):
1. Add the new variable as an "outvar" object in the "derivvarinit" function in dvardict.py.
2. Run dvardict.py (This will generate a new dictionary connecting the "verbose" variable names to the "rams-style" variable names)
3. Make a function in derivcalcs.py to calculate that derived variable (which should be in the units you specified in dvardict.py)
4. Add a code section in the "get_derivedvars" function to call the calculation function and append that variable to the "vardict" dictionary containing all other post-processed data.
4a. An example of this section is as follows:
    if "NEWVARRAMSNAME" in vardict.keys():
        newvard = get_newvar(vardict) #Depending on the variable you want to calculate, you may want to use the "rawfile" with sigma-z coordinates instead of the interpolated "vardict". For example, 
                                      #PWAT is  calculated using get_pwat(rawfile), since that uses the sigma-z coordinates for more accurate integration
        vardict["NEWVARRAMSNAME"].data = newvard
        print(f"Found variable {[vardict['NEWVARRAMSNAME'].varname, vardict['NEWVARRAMSNAME'].ramsname][rnameflag]}!")

>>>>>>> bf58ee02f99d24f00eba4863cdc1d5d72b7047d3
