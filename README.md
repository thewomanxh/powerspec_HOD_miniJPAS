These codes are to calculate Auto correlation power spectrum for miniJPAS based HOD.
Copyright (c) 2021, Xiu-hui Tan


## File Description ##

** <param_fac> stores parameters and arguments, which can be modified.

** <Func_Cosmo>: several functions for cosmological quantities, which includes **Cosmological Distance**, **Cosmological Time**, **Time vs Redshift convertion**, **comoving volume**, **dnde function of dark matter**, **optical depth**, et al.

** <powerspec_HOD> stores functions those to compute window function, number density, power spectrum from measurement and theory.
You can run MCMC with the functions in powerspec_HOD to find the HOD parameters.

** <read_miniJPAS_cl> is to build input data for powerspec_HOD.

** <read_miniJPAS_autodata> is to select data from the original one.

** <fit_HOD> stores functions those to fit hod parameters, you can get the results from <show_fitHOD>

** <show_fitHOD> is to calculate HOD parameters fitting results, which you may want to use them as an input lookup file in <hodpy>

** <divide_miniJPAS_field> is to build subregions of miniJPAS by dividing it into equal area field.


## Operating Instructions ##

1. check file path: /Data_Input/JPAS/minijpas_pdr201912_cat_wAbsMag_v1.csv should be exit, run <read_miniJPAS_data> to generate the band selection data in the /Data_Output/miniJPAS_selection/.

2. check file path: /Data_Input/Measurement_Cl_from_carlos/, run <read_miniJPAS_cl> to get Cl measurement Cl of miniJPAS from Carlos. The results in the /Data_Output/miniJPAS_cl/.

3. obtain Measurement and Theory Cl from <powerspec_HOD>, run MCMC to get bestfits of HOD parameters by <fit_Cl>. The results in the /Data_Output/Fit_res/.

4. build a file /Data_Output/Xiu_out/hod_3param.dat to store the output results from /Data_Output/Fit_res/HOD_*_3param.npy.

5. run <show_fitHOD> to get 3 HOD parameters fitting results, which the outputs are in the /Data_Output/HOD_params_res/. The results are used to be the insert parameters in /hodpy_miniJPAS/lookup/hod_params.dat. Find corresponding term and modify the values.






