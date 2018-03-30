Julia package for optimal spatial stratification

The package ospats+ contains 4 Julia script files:

1) “main”: the script for setting process parameters and calling the function ospats() or function ospall()

2) “readdata”: the script for reading the input file with grid data.

3) “ospats”: containing only the function ospats().

4) “ospall”: containing only the function ospall(). This function should be used if coarse-gridding is necessary, i.e. if the grid has too many points to enable processing of the pairwise distance matrix.

The functions ospats() and ospall() both lead to a stratification by the Ospats method. The difference is that opspats() applies the method to the entire grid, while ospall() applies the method to a systematic sample from the grid, followed by allocation of the remaining grid points to the strata as constructed from the sample.

The requested function is called by setting the parameter "in" (from sampling INterval) in the main script. If in=1 then ospats() is called, if “in” is set to a value larger than 1, then ospall() is called. For instance, if in=3, then every third grid point gets included in the sample.

For further explanation see the user manual at https://github.com//jjdegruijter/ospats+/Manual.