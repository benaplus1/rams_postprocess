
=======
# rams_postprocess
Post-Processing code for RAMS model output

Can interpolate variables to cartesian or pressure levels, or keep the data on a sigma-z grid
Interpolates vector quantities (like u, v, w) from grid cell edges of the Arikawa C-grid to the grid cell centers
Returns many derived variables from model output (such as pressure, temperature, dewpoint, particle number concentrations, particle diameters, buoyancy, and more!)
Parallelized using the builtin concurrent.futures package in Python
>>>>>>> bf58ee02f99d24f00eba4863cdc1d5d72b7047d3
