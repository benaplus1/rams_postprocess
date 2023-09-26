
=======
# rams_postprocess
Post-Processing code for RAMS model output

Can interpolate variables to cartesian or pressure levels, or keep the data on a sigma-z grid
Interpolates vector quantities (like u, v, w) from grid cell edges of the Arikawa C-grid to the grid cell centers
Returns many derived variables from model output (such as pressure, temperature, dewpoint, particle number concentrations, particle diameters, buoyancy, and more!)
Parallelized using the builtin concurrent.futures package in Python

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
