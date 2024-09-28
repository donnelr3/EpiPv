rm(list=setdiff(ls(), "filepath"))
##############################################################################################################
# This is the master script for the data underlying tables and figures in main text of                       #
# Donnelly Tankam Gilligan 2024 (henceforth MS). It has two main purposes. First it calls                    #
# estimate_virus_parameters_SPT and estimate_virus_parameters_PT - functions that estimate                   # 
# parameters from SPT and PT viral access period data as described in MS.                                    #
# Second it calls the calculate_epidemic_probability function which infers epidemic risk from PT             # 
# and SPT viral parameter estimates. The script exports summary data for the parameter estimates as          # 
# well the full Markov chains of posterior samples. The script reads in the Markov chains of                 # 
# posterior samples for the virus parameter estimates for the calculate_epidemic_probability                 # 
# function. The script exports summary data for epidemic risk as well as full chains for extinction risk     #
# depending on whether the virus arrived through infected insect imigration or infected plant propagation    #  
##############################################################################################################

setwd("C:/Users/ruair/Desktop/EpiPv")
#setwd("C:/Users/donnelr3/OneDrive - University Of Cambridge/Archive_ZG_WaveCassava/rFiles/rStanFiles/RPanalysis")

# At present source the following directory functions - But in final published forms this will be library('EpiPv')
source("solveInoculumStatesBP.r")
source("estimate_virus_parameters_PT.r")
source("estimate_virus_parameters_SPT.r")
source("calculate_epidemic_probability.r")

library(rstan) 
library(bridgesampling)
library(rstudioapi)
#library(nleqslv)
library(posterior)


options(mc.cores = parallel::detectCores()) 
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')


virus_type='PT'
run_PT_est_PUB=0
run_PT_est_SIM=0
run_SPT_est_PUB=0
calc_p_epi=1

nChains=4

################################################################################################################
################### DATA #######################################################################################
################################################################################################################
#### SAVE THIS AS A DATATABLE THAT LIVES IN THE PACKAGE AND IS LOADED HERE #####################################
# DUBERN RETENTION DATA (see Table 1 of corresponding Donnelly and GIlligan 2023 ms)
AAP_lens=c(2*60,3*60,3.5*60,4*60,4.5*60,5*60,6*60,8*60); # theres one less AP period hence taggging on a 0
LAP_lens=c(30,1*60,2*60,3*60,4*60,5*60,6*60,7*60,8*60);
IAP_lens=c(5,10,15,20,25,30,40,50,60); # MINUTES!

AAP_Reps=rep(30,length(AAP_lens)); # theres one less AP period hence taggging on a 0
LAP_Reps=c(12,12,22,22,34,34,34,34,34);
IAP_Reps=c(12,36,36,36,36,24,24,12,12);

AAP_Infs=c(0,0,9,13,11,16,21,16); # theres one less AP period hence taggging on a 0
LAP_Infs=c(0,0,0,0,3,11,19,21,23);
IAP_Infs=c(0,8,11,10,16,17,19,8,9);

lsEst_in=4*24*60;
# SEE ABOVE INSTRUCTION
# MARUTHI RETENTION DATA (see Table of )
AAP_lens_mg=c(10,30,60,60*4,60*24,60*48);  # HAVE ROUNDED FIRST ENTRY FROM 7.5!
LAP_lens_mg=c(60*48);
IAP_lens_mg=c(10,30,60,60*4,60*24,60*48); # MINUTES!

AAP_Reps_mg=c(25,25,25,15,20,15); # theres one less AP period hence taggging on a 0
LAP_Reps_mg=c(15);
IAP_Reps_mg=c(31,33,39,35,48,15);

AAP_Infs_mg=c(4,8,10,6,9,6); # theres one less AP period hence taggging on a 0
LAP_Infs_mg=c(0);
IAP_Infs_mg=c(6,7,8,13,29,6);

# SIMULATED PT DATASET
TPI_assay1_PT=as.matrix(read.table(('TPI_assay1_PT.dat')))
TPI_assay2_PT=as.matrix(read.table(('TPI_assay2_PT.dat')))
TPI_assay3_PT=as.matrix(read.table(('TPI_assay3_PT.dat')))

T_vecs_PT=(dget("assay123_dataList_PT.txt"))

AAP_lens_test_PT=unlist(T_vecs_PT[1]) # theres one less AP period hence taggging on a 0
LAP_lens_test_PT=unlist(T_vecs_PT[2])
IAP_lens_test_PT=unlist(T_vecs_PT[3])

AAP_Reps_test_PT=rep(dim(TPI_assay1_PT)[2],length(AAP_lens_test_PT)); # theres one less AP period hence taggging on a 0
LAP_Reps_test_PT=rep(dim(TPI_assay2_PT)[2],length(LAP_lens_test_PT));
IAP_Reps_test_PT=rep(dim(TPI_assay3_PT)[2],length(IAP_lens_test_PT));

AAP_Infs_test_PT=rowSums(TPI_assay1_PT); # theres one less AP period hence taggging on a 0
LAP_Infs_test_PT=rowSums(TPI_assay2_PT);
IAP_Infs_test_PT=rowSums(TPI_assay3_PT);

T_A_PT=unlist(T_vecs_PT[4])
T_L_PT=unlist(T_vecs_PT[5])
T_I_PT=unlist(T_vecs_PT[6])
numFly_PT=unlist(T_vecs_PT[7])
################################################################################################################
################### ESTIMATES ##################################################################################
################################################################################################################
#########   pass  table to estimate_virus_parameters together with D_WF0 and D_bgLens for parameter ESTIMATES
# DUBERN PT #
if (run_PT_est_PUB) {
  EVPT_pub=estimate_virus_parameters_PT(rbind(AAP_lens,AAP_Reps,AAP_Infs),
                                        rbind(LAP_lens,LAP_Reps,LAP_Infs),
                                        rbind(IAP_lens,IAP_Reps,IAP_Infs),
                                        rbind(c(0,0,24*3.5*60),c(5*60,0,24*3.5*60),c(5*60,6*60,0)),
                                        10,lsEst_in,'CMB')
  # nb estimate_virus_parameters_PT by default exports a summary table of parameter estimates
  saveRDS(EVPT_pub, file="EVPT_pub") # in addition exporting chains
}else{
  # here's one i made earlier, estimates for published PT virus
  EVPT_pub <- readRDS("EVPT_pub")
}
# SIMULATED PT #
if (run_PT_est_SIM) {
  EVPT_sim=estimate_virus_parameters_PT(rbind(AAP_lens_test_PT,AAP_Reps_test_PT,AAP_Infs_test_PT),
                                        rbind(LAP_lens_test_PT,LAP_Reps_test_PT,LAP_Infs_test_PT),
                                        rbind(IAP_lens_test_PT,IAP_Reps_test_PT,IAP_Infs_test_PT),
                                        rbind(c(0,T_L_PT,T_I_PT),c(T_A_PT,0,T_I_PT),c(T_A_PT,T_L_PT,0)),
                                        numFly_PT,lsEst_in,'SIM')
  # nb estimate_virus_parameters_PT by default exports a summary table of parameter estimates
  saveRDS(EVPT_sim, file="EVPT_sim") # in addition exporting chains
}else{
  # here's one i made earlier, estimates for simulated PT virus
  EVPT_sim <- readRDS("EVPT_sim")
}
# MARUTHI SPT #
if (run_SPT_est_PUB) {
  EVSPT_pub=estimate_virus_parameters_SPT(rbind(AAP_lens_mg,AAP_Reps_mg,AAP_Infs_mg),
                                          rbind(IAP_lens_mg,IAP_Reps_mg,IAP_Infs_mg),
                                          rbind(rbind(c(0,60*48),c(60*48,0))),
                                          23,lsEst_in,'CBSI')
  # nb estimate_virus_parameters_PT by default exports a summary table of parameter estimates
  saveRDS(EVSPT_pub, file="EVSPT_pub") # in addition exporting chains
}else{
  # here's one i made earlier, estimates for published SPT virus
  EVSPT_pub <- readRDS("EVSPT_pub")
}


if (virus_type=='PT') {
  target=EVPT_pub; # setting the chains in question for subsequent re-use
  numParams=5;
}else if(virus_type=='SPT'){
  target=EVSPT_pub; # setting the chains in question for subsequent re-use
  numParams=3;
}else{
  print('did not set for valid plant virus type')
}

target_B=array(rep(0, (dim(target)[1]/nChains)), dim=c((dim(target)[1]/nChains), nChains, numParams))
for (gg in 1:nChains){
  target_B[,gg,]=target[(((gg-1)*(dim(target)[1]/nChains))+1):((gg)*(dim(target)[1]/nChains)),]
  colnames(target_B[,gg,])=colnames(target)
}

# FOR INFO
tmpb=as_draws(target_B)
sdraws=summarise_draws(tmpb)
# can check against stan parameter summaries
# may need to convert from mns to hrs
################################################################################################################
################################################################################################################



################################################################################################################
################### INFERENCES #################################################################################
################################################################################################################
# P EPIDEMIC inferences
# accessing the chains from estimate_virus_parameters 
# VIRUS PARAMETERS
al_fits=target[,which(colnames(target)=='al[1]')]*60*24 # convert from per min to per day
be_fits=target[,which(colnames(target)=='be[1]')]*60*24
if (virus_type=='PT') {
  mu_fits=target[,which(colnames(target)=='mu[1]')]*60*24
  # HOWEVER - for PT only mu represents rate of loss of insect infectiousness which may be through mortality of clearance
  # user-defined lsEst_in provides estimate of laboratory insect survival. target for PT contains mu_cl[1]: separated out virus clearance
  # for CMB Dubern dataset the median mu_cl[1] is very close to zero representing a probably absence of B. tabaci clearancce of CMB
  # Accordingly:
  mu_fits=target[,which(colnames(target)=='mu[1]')]*0
  lat_fits=target[,which(colnames(target)=='lat[1]')]*60*24 # note that latent progression not currently used in field model as relatively short
  }else{
  mu_fits=target[,which(colnames(target)=='mu[1]')]*60*24
}
# sample 1000 sets of parameters for calculating epidemic probability
runs=10    # length(al_fits)
# LOCAL PARAMETERS
# set the local parameters (per day)
thet_external <- 0.45 # dispersal 
r_external  <- 1/21 # roguing 
bf_external <- 1/21  # vector mortality 
h_external <- 1/365  # harvesting rate 
nu_pl_external <- 1/14  # plant latent period 
localParams=c(thet_external, r_external, h_external, bf_external, nu_pl_external)
# generate p Epdemic results for various insect burdens
numInsects_vec_cbsi=c(4,6,8)
numInsects_vec_cmb=c(1,2,3)
if (virus_type=='PT') {
  numInsects_vec=numInsects_vec_cmb; 
}else if(virus_type=='SPT'){
  numInsects_vec=numInsects_vec_cbsi; 
}else{
  print('did not set for valid plant virus type')
}

if (calc_p_epi) {
  data_table_Pl=matrix(0, length(numInsects_vec), 3)
  data_table_Ins=matrix(0, length(numInsects_vec), 3)
  pepi_frPlant=matrix(0, length(numInsects_vec), runs) 
  pepi_frInsect=matrix(0, length(numInsects_vec), runs)
  for (ii in 1:3){ #  <---------- set the insect burden for calculating epidemic probability
    numInsects=numInsects_vec[ii]
    print('numInsects')
    print(numInsects)
    numVars <- ((numInsects+1)*3)-1 
    result_vec <- matrix(0, 2, runs)  
    interval_ind=10
    # loop through MCMC chains to sample posterior distributions of virus parameters
    for (ppp in 1:runs) {  #
      print(ppp)
      al_estim=al_fits[length(al_fits)+1-ppp]
      be_estim=be_fits[length(al_fits)+1-ppp]
      mu_estim=mu_fits[length(al_fits)+1-ppp]
      #mu_estim=lat_fits[ppp] vector latency not currently used in field model as vector latency is relatively short hence not significant
      virusParams=c(al_estim,be_estim,mu_estim)
      
      # calculate_epidemic_probability_vBP returns vector of epidemic probabilities for differnt inoculum states
      qm_out=calculate_epidemic_probability(numInsects,interval_ind,localParams,virusParams)
      result_vec[1,(ppp)]=qm_out[1] 
      result_vec[2,(ppp)]=qm_out[(numVars-(numInsects-1))] 
    }
    pepi_frPlant[ii,]=result_vec[1,] # storing extinction probabilities from plant inoculum
    pepi_frInsect[ii,]=result_vec[2,] # and from vector inoculum
    
    # for a given level of insect summarise_draws for credible intervals
    target_A=array(rep(0, (runs/nChains)), dim=c((runs/nChains), nChains, 2))
    for (gg in 1:nChains)
      target_A[,gg,]=t(result_vec[,(((gg-1)*(runs/nChains))+1):((gg)*(runs/nChains))])
    
    ts=summarize_draws(as_draws(target_A))
    data_table_Pl[ii,]=c(ts$median[1],ts$q5[1],ts$q95[1])
    data_table_Ins[ii,]=c(ts$median[2],ts$q5[2],ts$q95[2])
  }
  #### EXPORT TABLE
  # FIRST EXPORT ENTIRE epidemic probability CHAINS
  ###################################################################
  row.names(pepi_frPlant)=numInsects_vec
  write.table(pepi_frPlant, paste('path_',virus_type,'_chains_byF_frPlant.dat'), append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
  row.names(pepi_frInsect)=numInsects_vec
  write.table(pepi_frInsect, paste('path_',virus_type,'_chains_byF_frInsect.dat'), append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
  ###################################################################
  
  
  #### export summary tables for p epidemic INFERENCES ####
  #data_table_Pl=data_table_Pl*(data_table_Pl>=0)
  #data_table_Pl=data_table_Pl*(data_table_Pl<=1)+1*(data_table_Pl>1)
  row.names(data_table_Pl)=numInsects_vec
  colnames(data_table_Pl)=c('50','2.5','97.5')
  write.table(data_table_Pl, paste('path_',virus_type,'_summary_byF_frPlant.dat'), append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
  
  #data_table_Ins=data_table_Ins*(data_table_Ins>=0)
  #data_table_Ins=data_table_Ins*(data_table_Ins<=1)+1*(data_table_Ins>1)
  row.names(data_table_Ins)=numInsects_vec
  colnames(data_table_Ins)=c('50','2.5','97.5')
  write.table(data_table_Ins, paste('path_',virus_type,'_summary_byF_frInsect.dat'), append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
}else{
  # here's one i made earlier, P EPIDEMICS chains
  old_chaindata_table_Pl=as.matrix(read.table((paste('path_',virus_type,'_chains_byF_frPlant.dat'))))
  old_chaindata_table_Ins=as.matrix(read.table((paste('path_',virus_type,'_chains_byF_frInsect.dat'))))
  # here's one i made earlier, P EPIDEMICS summaries
  old_data_table_Pl=as.matrix(read.table(('pathIntervals_byNumIns_frPlant.dat')))
  old_data_table_Ins=as.matrix(read.table(('pathIntervals_byNumIns_frIns.dat')))
}

################################################################################################################
################################################################################################################