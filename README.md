These codes are to calculate Auto correlation power spectrum for miniJPAS based HOD.


** <powerspec_HOD> stores functions those to compute window function, number density, power spectrum from measurement and theory.
You can run MCMC with the functions in powerspec_HOD to find the HOD parameters.

** <param_fac> stores parameters, which can be modified.

** <read_miniJPAS_cl> is to build input data for powerspec_HOD.

** <read_miniJPAS_autodata> is to select data from the original one.

** <fit_HOD> stores functions those to fit hod parameters, you can get the results from <show_fitHOD>

** <show_fitHOD> is to calculate HOD parameters fitting results, which you may want to use them as an input lookup file in <hodpy>

** <divide_miniJPAS_field> is to build subregions of miniJPAS by dividing it into equal area field.


