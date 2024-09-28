# EpiPv
R package functions for estimating plant virus transmission parameters and epidemic risk

After agreeing on a final version of a corresponding manuscript with a prospective academic journal updated versions of the function files found here will be released as a R package.

The main functions are:

The {\it estimate\_virus\_parameters\_PT} function receives AP data for a given vector-PT virus-plant combination as user input together with assay configuration (i.e., AP feeding durations $T_A$ $T_L$ $T_I$ and the number of insect vectors used $X_0$) and returns posterior parameter distributions for the transmission rates $\mu$, $\alpha$, $\beta$, and $\gamma$ (see Appendix S1 for details).

The {\it estimate\_virus\_parameters\_SPT} function receives AP data for a given vector-SPT virus-plant combination as user input together with assay configuration (i.e., AP feeding durations $T_A$ $T_I$ and the number of insect vectors used $X_0$) and returns posterior parameter distributions for the rates transmission rates $\mu$, $\alpha$, and $\beta$ (see Appendix S2 for details).

The {\it calculate\_epidemic\_probability} function receives event rate parameters for virus transmission ($\mu$, $\alpha$, $\beta$) as well as local parameters $F$, $\theta$, $r$, $h$ and $b_f$ as user input, and returns epidemic probability for different types of inoculum state (see Appendix S3 for details).

The {\it AP\_data\_simulator} function receives event rate parameters for virus transmission ($\mu$, $\alpha$, $\beta$) as well as assay feeding durations and returns simulated access period data (see Appendix S4 for details of the statistical simulation process).

Additional functions and scripts are:
- APdata_simulation_script (it sets the assay structure and calls the AP\_data\_simulator function)
- inoc_durtn_calculator (it calculates the effective duration an insect vector has for inoculating a test plant)
- model_AP_PT.stan (it encodes parameter estimation though probability model of the AP assay for a PT virus)
- model_AP_SPT.stan (it encodes parameter estimation though probability model of the AP assay for an SPT virus)
- ppath_profilr_masterscript_BP (it sets the assay data from which to estimate parameters; it calls estimate\_virus\_parameters\_PT, estimate\_virus\_parameters\_SPT and calculate\_epidemic\_probability functions)
- solveInoculumStatesBP (this function receives virus parameters and local parameters and returns two matrices that are used to infer epidemic risk)
  
Additional data files are in the data_files folder - this includes Markov chain results from the main paper for a real PT virus, a simulated PT virus and a real SPT virus that can be called in ppath_profilr_masterscript_BP to infer epidemic risk
