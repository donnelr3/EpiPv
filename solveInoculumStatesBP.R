solveInoculumStatesBP <- function(numVars, ext_pams, est_pams,interv) {
  
  # this code represents the system for extinction probability of invading insect-borne plant pathogens
  # this function returns two rbinded arrays ttilde (population transitions) and f (fertility) as outlined in
  # Donnelly, Tankam and Gilligan 2024 (henceforth MS)
  # the challenge here is to do this flexibly as different levels of the constant insect burden per plant
  # will produce different dimensions of ttilde and f
  
  # MINI FUNCTION
  rateToProb <- function(eventRate,intervalL) {
    eventProb=1-exp(-eventRate*intervalL)
    return(eventProb)
  }
  
  # incoming rates must be in a per day format
  #print(ext_pams)
  thet_det <- ext_pams[1] # dispersal
  r_det <- ext_pams[2] # roguing
  h_det <- ext_pams[3] # harvest
  vm_det <- ext_pams[4] # vector mortality
  mu_pl_det <- ext_pams[5]  # plant latent
  #print(est_pams)
  al_est <- est_pams[1] # acquire
  be_est <- est_pams[2] # inoculate
  mu_est <- est_pams[3]  # retain
  #mu_est <- est_pams[4]  # vector latent
  
  numIns=((numVars+1)/3)-1;
  
  if (!numIns > 0) {
    stop(paste("Error user inputs 0 insect vectors per plant:", numIns))
  }else{
    paste("User inputs insect vectors per plant:", numIns)
  }
  
  A_mat <- matrix(0, numVars, numVars)
  B_vec <- rep(0, numVars)
  F_mat <- matrix(0, numVars, numVars)
  
  # useful to assign index numbers to inoculum states as the challenge here is dimension flexibility 
  # (as different levels of the constant insect burden per plant will produce different arrray dimensions which we are then populating)
  indI_0_wf <- 1 + 0 * (numIns + 1) + 0
  indI_F_wf <- 1 + 0 * (numIns + 1) + numIns
  indE_0_wf <- 1 + 1 * (numIns + 1) + 0
  indE_F_wf <- 1 + 1 * (numIns + 1) + numIns
  indS_F_wf <- 2 * (numIns + 1) + numIns
  
  # indices for intermediate infected insect values per inoculum state
  if (numIns > 1) {
    indI_wf <- 1 + 0 * (numIns + 1) + (1:(numIns - 1))
    indE_wf <- 1 + 1 * (numIns + 1) + (1:(numIns - 1))
    indS_wf <- 2 * (numIns + 1) + (1:(numIns - 1))
  }
  
  #  block I[1], rates applying to infected plant inoculum state with no infected insects
  if (numIns == 1) {
    A_mat[indI_F_wf,indI_0_wf ] <- A_mat[ indI_F_wf,indI_0_wf] + rateToProb((numIns * al_est ),interv)
  } else {
    A_mat[indI_wf[1],indI_0_wf ] <- A_mat[indI_wf[1],indI_0_wf ] + rateToProb((numIns * al_est ),interv)
  }
  B_vec[indI_0_wf] <- B_vec[indI_0_wf] + rateToProb(((h_det + r_det) ),interv)
  
  # block E[1], rates applying to exposed plant inoculum state with no infected insects
  A_mat[indI_0_wf,indE_0_wf ] <- A_mat[indI_0_wf,indE_0_wf ] + rateToProb((mu_pl_det ),interv)
  B_vec[indE_0_wf] <- B_vec[indE_0_wf] + rateToProb((h_det ),interv)

  if (numIns > 1) {
    
    #  block I[1+kkk], rates applying to infected plant inoculum state with kkk infected insects
    for (kkk in 1:(numIns - 1)) {
      if (kkk < (numIns - 1)) {
        A_mat[indI_wf[kkk + 1],indI_wf[kkk] ] <- A_mat[indI_wf[kkk + 1],indI_wf[kkk] ] + rateToProb(((numIns - kkk) * al_est ),interv)
      } else {
        A_mat[indI_F_wf,indI_wf[kkk] ] <- A_mat[indI_F_wf,indI_wf[kkk] ] + rateToProb((al_est ),interv)
      }
      if (kkk > 1) {
        A_mat[indI_wf[kkk - 1],indI_wf[kkk] ] <- A_mat[indI_wf[kkk - 1],indI_wf[kkk] ] + rateToProb((kkk * (vm_det + mu_est) ),interv)
      } else {
        A_mat[indI_0_wf,indI_wf[1] ] <- A_mat[indI_0_wf,indI_wf[1] ] + rateToProb((1 * (vm_det + mu_est) ),interv)
      }
      A_mat[indS_wf[kkk],indI_wf[kkk] ] <- A_mat[indS_wf[kkk],indI_wf[kkk] ] + rateToProb(((h_det + r_det) ),interv)
      
      
      #  block E[1+kkk], rates applying to exposed plant inoculum state with kkk infected insects 
      A_mat[indI_wf[kkk],indE_wf[kkk] ] <- A_mat[indI_wf[kkk],indE_wf[kkk] ] + rateToProb((mu_pl_det ),interv)
      if (kkk > 1) {
        A_mat[indE_wf[kkk - 1],indE_wf[kkk] ] <- A_mat[indE_wf[kkk - 1],indE_wf[kkk] ] + rateToProb((kkk * (vm_det + mu_est) ),interv)
      } else {
        A_mat[indE_0_wf,indE_wf[1] ] <- A_mat[indE_0_wf,indE_wf[1] ] + rateToProb((1 * (vm_det + mu_est) ),interv)
      }
      A_mat[indS_wf[kkk],indE_wf[kkk] ] <- A_mat[indS_wf[kkk],indE_wf[kkk] ] + rateToProb((h_det ),interv)
      
      
      # block S[1+kkk], rates applying to susceptible plant inoculum state with kkk infected insects 
      A_mat[indE_wf[kkk],indS_wf[kkk] ] <- A_mat[indE_wf[kkk],indS_wf[kkk] ] + rateToProb((kkk * be_est ),interv)
      if (kkk > 1) {
        A_mat[indS_wf[kkk - 1],indS_wf[kkk] ] <- A_mat[indS_wf[kkk - 1],indS_wf[kkk] ] + rateToProb((kkk * (vm_det + mu_est) ),interv)
      } else {
        B_vec[indS_wf[1]] <- B_vec[indS_wf[1]] + rateToProb((1 * (vm_det + mu_est)),interv)
      }
    }
  }
  
  
  # block I[1+numIns], rates applying to infected plant inoculum state with all numIns insects infected
  if (numIns == 1) {
    A_mat[indI_0_wf,indI_F_wf ] <- A_mat[indI_0_wf,indI_F_wf ] + rateToProb((numIns * (vm_det + mu_est) ),interv)
  } else {
    A_mat[indI_wf[numIns - 1],indI_F_wf ] <- A_mat[indI_wf[numIns - 1],indI_F_wf ] + rateToProb((numIns * (vm_det + mu_est) ),interv)
  }
  A_mat[indS_F_wf,indI_F_wf ] <- A_mat[indS_F_wf,indI_F_wf ] + rateToProb(((h_det + r_det) ),interv)
  
  
  # block E[1+numIns], rates applying to exposed plant inoculum state with all numIns insects infected
  A_mat[indI_F_wf,indE_F_wf ] <- A_mat[indI_F_wf,indE_F_wf ] + rateToProb((mu_pl_det ),interv)
  A_mat[indS_F_wf,indE_F_wf ] <- A_mat[indS_F_wf,indE_F_wf ] + rateToProb((h_det ),interv)
  
  
  # block S[1+numIns], rates applying to susceptible plant inoculum state with all numIns insects infected
  A_mat[indE_F_wf,indS_F_wf ] <- A_mat[indE_F_wf,indS_F_wf ] + rateToProb((numIns * be_est),interv)
  
  
  # finally ensuring entries in the right place: pay particular intention to case of no intermediate insects (only 0 or all infected i.e. when numIns=1)
  if (numIns==1) {
    A_mat[indE_0_wf,indE_F_wf ] <- A_mat[indE_0_wf,indE_F_wf ] + rateToProb((numIns * (vm_det + mu_est) ),interv)
    B_vec[indS_F_wf] <- B_vec[indS_F_wf] + rateToProb((1 * (vm_det + mu_est) ),interv)
  }else{
    A_mat[indE_wf[numIns-1],indE_F_wf ] <- A_mat[indE_wf[numIns-1],indE_F_wf ] + rateToProb((numIns * (vm_det + mu_est) ),interv)
    A_mat[indS_wf[numIns-1],indS_F_wf ] <- A_mat[indS_wf[numIns-1],indS_F_wf ] + rateToProb((numIns * (vm_det + mu_est) ),interv)
  }
  

  
  # finally pay special attention to infected insect dispersal - in addition to population transitions that result (A_mat) this is
  # the only event that leads to 'fertility' (since at invasion they will land on healthy plants and see new inocula)
  if (numIns == 1) {
    A_mat[indI_0_wf,indI_F_wf ] <- A_mat[indI_0_wf,indI_F_wf ] + rateToProb(((thet_det)),interv) # LOSE INSECT P ext IYj->P ext IYj-1
    A_mat[indE_0_wf,indE_F_wf ] <- A_mat[indE_0_wf,indE_F_wf ] + rateToProb(((thet_det)),interv) # LOSE INSECT P ext EYj->P ext EYj-1
    
    B_vec[indS_F_wf] <- B_vec[indS_F_wf] + rateToProb(((thet_det)),interv)

    F_mat[indS_F_wf,indI_F_wf ] <- F_mat[indS_F_wf,indI_F_wf ] + rateToProb(((thet_det)),interv) # LOSE INSECT P ext IYj->P ext IYj-1
    F_mat[indS_F_wf,indE_F_wf ] <- F_mat[indS_F_wf,indE_F_wf ] + rateToProb(((thet_det)),interv) # LOSE INSECT P ext EYj->P ext EYj-1
    F_mat[indS_F_wf,indS_F_wf ] <- F_mat[indS_F_wf,indS_F_wf ] + rateToProb(((thet_det)),interv) # LOSE INSECT P ext SYj->P ext SYj-1
    
  } else {
    A_mat[indI_0_wf,indI_wf[1] ] <- A_mat[indI_0_wf,indI_wf[1] ] + rateToProb(((thet_det)),interv) # LOSE INSECT P ext IYj->P ext IYj-1
    A_mat[indE_0_wf,indE_wf[1] ] <- A_mat[indE_0_wf,indE_wf[1] ] + rateToProb(((thet_det)),interv) # LOSE INSECT P ext EYj->P ext EYj-1
    
    B_vec[indS_wf[1]] <- B_vec[indS_wf[1]] + rateToProb(((thet_det)),interv)

    F_mat[indS_wf[1],indI_wf[1] ] <- F_mat[indS_wf[1],indI_wf[1] ] + rateToProb(((thet_det)),interv) # LOSE INSECT P ext IYj->P ext IYj-1
    F_mat[indS_wf[1],indE_wf[1] ] <- F_mat[indS_wf[1],indE_wf[1] ] + rateToProb(((thet_det)),interv) # LOSE INSECT P ext EYj->P ext EYj-1
    F_mat[indS_wf[1],indS_wf[1] ] <-  F_mat[indS_wf[1],indS_wf[1] ] + rateToProb(((thet_det)),interv) # LOSE INSECT P ext SYj->P ext SYj-1
    
    for (ww in 2:(numIns-1)) {  
      A_mat[indI_wf[ww-1],indI_wf[ww] ] <- A_mat[indI_wf[ww-1],indI_wf[ww] ] + rateToProb(ww*((thet_det)),interv) # LOSE INSECT P ext IYj->P ext IYj-1
      A_mat[indE_wf[ww-1],indE_wf[ww] ] <- A_mat[indE_wf[ww-1],indE_wf[ww] ] + rateToProb(ww*((thet_det)),interv) # LOSE INSECT P ext EYj->P ext EYj-1
      A_mat[indS_wf[ww-1],indS_wf[ww] ] <- A_mat[indS_wf[ww-1],indS_wf[ww] ] + rateToProb(ww*((thet_det)),interv) # LOSE INSECT P ext SYj->P ext SYj-1
      
      F_mat[indS_wf[1],indI_wf[ww] ] <- F_mat[indS_wf[1],indI_wf[ww] ] + rateToProb(ww*((thet_det)),interv) # LOSE INSECT P ext IYj->P ext IYj-1
      F_mat[indS_wf[1],indE_wf[ww] ] <- F_mat[indS_wf[1],indE_wf[ww] ] + rateToProb(ww*((thet_det)),interv) # LOSE INSECT P ext EYj->P ext EYj-1
      F_mat[indS_wf[1],indS_wf[ww] ] <- F_mat[indS_wf[1],indS_wf[ww] ] + rateToProb(ww*((thet_det)),interv) # LOSE INSECT P ext SYj->P ext SYj-1
    }
    A_mat[indI_wf[numIns-1],indI_F_wf ] <- A_mat[indI_wf[numIns-1],indI_F_wf ] + rateToProb(numIns*((thet_det)),interv) # LOSE INSECT P ext IYj->P ext IYj-1
    A_mat[indE_wf[numIns-1],indE_F_wf ] <- A_mat[indE_wf[numIns-1],indE_F_wf ] + rateToProb(numIns*((thet_det)),interv) # LOSE INSECT P ext EIYj->P ext EYj-1
    A_mat[indS_wf[numIns-1],indS_F_wf ] <- A_mat[indS_wf[numIns-1],indS_F_wf ] + rateToProb(numIns*((thet_det)),interv) # LOSE INSECT P ext IYj->P ext IYj-1
    
    F_mat[indS_wf[1],indI_F_wf ] <- F_mat[indS_wf[1],indI_F_wf ] + rateToProb(numIns*((thet_det)),interv) # LOSE INSECT P ext IYj->P ext IYj-1
    F_mat[indS_wf[1],indE_F_wf ] <- F_mat[indS_wf[1],indE_F_wf ] + rateToProb(numIns*((thet_det)),interv) # LOSE INSECT P ext EIYj->P ext EYj-1
    F_mat[indS_wf[1],indS_F_wf ] <- F_mat[indS_wf[1],indS_F_wf ] + rateToProb(numIns*((thet_det)),interv) # LOSE INSECT P ext IYj->P ext IYj-1
  }
  
  
  Apr_mat=rbind(A_mat,B_vec)
  
  
  # small interval length means that the sum of probabilities of events occurring do not exceed 1 - so we also 
  # populate the nothing happens entry with 1- all the other rates
  for (jj in 1:numVars) {
    Apr_mat[jj, jj] <- 1 - sum(Apr_mat[, jj])
  }
  
  return(rbind(Apr_mat,F_mat))
}