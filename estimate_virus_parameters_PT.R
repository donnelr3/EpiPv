# This is the main function call that in turn runs the stan file for PT virus parameter estimation.
# Though the function receives model event parameters and underlying assay structure as arguments the settings for rstan to use are specified here
# In addition, the parameter estimates as percentile summaries and the full MArkov chains are all exported from this function  
# to the main home directory 

estimate_virus_parameters_PT <- function(assay1,assay2,assay3,offdiag_array,IVD,lsEst,virus){
  
  numWarm=500
  numIter=1000
  adaptVal=1-(10^-1)
  numChains=4
  treeDepth=15

  AAP_lens=assay1[1,];
  LAP_lens=assay2[1,];
  IAP_lens=assay3[1,];
  
  AAP_Reps=assay1[2,];
  LAP_Reps=assay2[2,];
  IAP_Reps=assay3[2,];
  
  AAP_Infs=assay1[3,];
  LAP_Infs=assay2[3,];
  IAP_Infs=assay3[3,];
  
  ##### check diags of offdiag_array #####
  if (sum(diag(offdiag_array)>0)!=0){
    print('input error!');
  }
  
  # assembling all the data inputs for the stan estimation
  dat3= list(D_NumGrps=c(length(AAP_Reps),length(LAP_Reps),length(IAP_Reps)),#length(AAP_lens),
             D_Wf0=IVD, #,
             D_Ls=lsEst,
             D_LensA=AAP_lens,
             D_LensL=LAP_lens,
             D_LensI=IAP_lens,
             D_RepsA=AAP_Reps,
             D_RepsL=LAP_Reps,
             D_RepsI=IAP_Reps,
             D_InfsA=AAP_Infs,
             D_InfsL=LAP_Infs,
             D_InfsI=IAP_Infs,
             D_bgLens=offdiag_array)
  
  #### it may help to run with initial values but one can use stan defaults which may generate warnings for initial chain values
  # but it is expected to settle down after intial jumps - otherwsise one can specify initials as per below example
  #initVal1=0.001;
  #initVal2=0.01; 
  #initf1_dat3 <- function() list(mu=array(initVal1,1),
   #                              al=array(initVal1,1),
  #                               mu_lat=array(initVal2,1),
   #                              mu_be=array(30,1))
  
  #### run stan on the data as a model OF AP assays for PT data ################
  fitB1 = stan(file = "model_AP_PT.stan",
               data = dat3,
               #init=initf1_dat3,
               iter = numIter, #numIter,1, #
               warmup = numWarm, #numWarm,0, #
               chains=numChains, #1,
               control = list(max_treedepth = treeDepth,adapt_delta=adaptVal,stepsize=0.01))
  stanOutB1=summary(fitB1)
  stanDetailB1=stanOutB1$summary
  
  #### extract summary parameter estimates ####
  al_out=summary(fitB1,'al[1]')$summary[c(1,4,6,8)]
  be_out=summary(fitB1,'be[1]')$summary[c(1,4,6,8)]
  lat_out=summary(fitB1,'lat[1]')$summary[c(1,4,6,8)]
  muloss_out=summary(fitB1,'mu[1]')$summary[c(1,4,6,8)]
  # note that muloss_out estimate actually represents total rate of loss of insect infectiousness
  # but the data input to stan takes an estimate of insect survival in the lab
  # this is used to separate out true virus clearance from mortality: that is mu_out
  mu_out=summary(fitB1,'mu_cl[1]')$summary[c(1,4,6,8)]



  #### export summary table for parameter estimates - switch to hourly units ####
  cmb_estimates=rbind(c(al_out*60),c(be_out*60),c(mu_out*60),c(lat_out*60),c(muloss_out*60))
  row.names(cmb_estimates)=c('al','be','mu','lat','mu_cl')
  colnames(cmb_estimates)=c('mean','2.5','50','97.5')
  write.table(cmb_estimates, paste(virus,'_PT_path_estimates.dat'), append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
  
  #### return chains for main parameters ####
  params<- as.matrix(fitB1, pars = c("al[1]", "be[1]", "mu_cl[1]","lat[1]","mu[1]")) 
  return(params);
  
}