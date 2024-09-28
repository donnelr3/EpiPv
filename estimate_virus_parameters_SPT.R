# This is the main function call that in turn runs the stan file for SPT virus parameter estimation.
# Though the function receives model event parameters and underlying assay structure as arguments the settings for rstan to use are specified here
# In addition, the parameter estimates as percentile summaries and the full MArkov chains are all exported from this function  
# to the main home directory 


estimate_virus_parameters_SPT <- function(assay1,assay2,offdiag_array,IVD,lsEst,virus){
  
  numWarm=4500
  numIter=6000
  adaptVal=1-(10^-1)
  numChains=4
  treeDepth=15

  AAP_lens=assay1[1,];
  IAP_lens=assay2[1,];
  
  AAP_Reps=assay1[2,];
  IAP_Reps=assay2[2,];
  
  AAP_Infs=assay1[3,];
  IAP_Infs=assay2[3,];
  
  ##### check diags of offdiag_array #####
  if (sum(diag(offdiag_array)>0)!=0){
    print('input error!');
  }
  
  # assembling all the data inputs for the stan estimation
  dat4b= list(D_NumGrps=c(length(AAP_lens),length(IAP_lens)),#length(AAP_lens),length(LAP_Reps_mg)
              D_Wf0=IVD, #,
              D_Ls=lsEst,
              D_LensA=AAP_lens,
              D_LensI=IAP_lens,
              D_RepsA=AAP_Reps,
              D_RepsI=IAP_Reps,
              D_InfsA=AAP_Infs,
              D_InfsI=IAP_Infs,
              D_bgLens=offdiag_array)
  
  #### it may help to run with initial values but one can use stan defaults which may generate warnings for initial chain values
  # but it is expected to settle down after intial jumps - otherwsise one can specify initials as per below example
  # initVal1=0.01;
  # initVal2=0.1;
  # initf1_dat4 <- function() list(c2=array(initVal1,1),
                                 #c3=array(initVal1,1),
                                #c1=array(initVal1,1))
  

  #### run stan on the data as a model OF AP assays for SPT data ################
  fitB2 = stan(file = "model_AP_SPT.stan",
               data = dat4b,
               #init=initf1_dat4,
               iter = numIter, #numIter,1, #
               warmup = numWarm, #numWarm,0, #
               chains=numChains, #1,
               control = list(max_treedepth = treeDepth,adapt_delta=adaptVal,stepsize=0.01))
  stanOutB2=summary(fitB2)
  stanDetailB2=stanOutB2$summary
  
  #### extract summary parameter estimates ####
  c2_out=summary(fitB2,'c2[1]')$summary[c(1,4,6,8)]
  c3_out=summary(fitB2,'c3[1]')$summary[c(1,4,6,8)]
  c1_out=summary(fitB2,'c1[1]')$summary[c(1,4,6,8)]
  al_out=summary(fitB2,'al[1]')$summary[c(1,4,6,8)]
  be_out=summary(fitB2,'be[1]')$summary[c(1,4,6,8)]
  mu_out=summary(fitB2,'mu[1]')$summary[c(1,4,6,8)]
  albe_out=summary(fitB2,'albe[1]')$summary[c(1,4,6,8)]
  ### for SPT the al_out,  be_out and  mu_out estimates are actually approximations (they are derived from the compunt parameters)
  # one can attempt a second round of improving approximations (below) as per appendix S2 main text
  al_p_out=summary(fitB2,'al_p[1]')$summary[c(1,4,6,8)]
  be_p_out=summary(fitB2,'be_p[1]')$summary[c(1,4,6,8)]
  mu_p_out=summary(fitB2,'mu_p[1]')$summary[c(1,4,6,8)]
 

  #### export summary table for parameter estimates - switch to hourly units ####
  path_estimates=rbind(c(al_out*60),c(be_out*60),c(mu_out*60),c(al_p_out*60),c(be_p_out*60),c(mu_p_out*60))
  rownames(path_estimates)=c('al','be','mu','alp','bep','mup')
  colnames(path_estimates)=c('mean','2.5','50','97.5')
  write.table(path_estimates, paste(virus,'_SPT_path_estimates.dat'), append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)

  #### return chains for main parameters ####
  params<- as.matrix(fitB2, pars = c("al[1]", "be[1]", "mu[1]")) #head(params)
  return(params);
  
}