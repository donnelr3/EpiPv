// Description: Probability model of SPT access period assay used in Bayesian analysis, see Donnelly Tankam Gilligan 2024 appendix S2 
// data : access period data NOTE that we assume there are only acquisition access and inoculation access components for SPT virus
//
// infer alpha beta mu from data 
// note convention throughout is that variables with D_ in front all involve data coming as input

functions {
  
  
}

// DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA 
data {                                    // Data block
int<lower=0> D_NumGrps[2];  // number of groups to be studied in each of AAP IAP variants
int<lower=1> D_Wf0;        // ultimate number of whitefly ON TEST PLANTS
real<lower=1> D_Ls;       //expected insect lifespan in minutes UNUSED for SPT
int<lower=0> D_LensA[D_NumGrps[1]];        // Number of reps for each of AAP LAP IAP variants across studied groups
int<lower=0> D_LensI[D_NumGrps[2]];
int<lower=0> D_RepsA[D_NumGrps[1]];        // Number of reps for each of AAP LAP IAP variants across studied groups
int<lower=0> D_RepsI[D_NumGrps[2]]; 
int<lower=0> D_InfsA[D_NumGrps[1]];                        // Response variable (if transm)
int<lower=0> D_InfsI[D_NumGrps[2]];  
int<lower=0> D_bgLens[2,2];
}

// TRANSFORMED DATA TRANSFORMED DATA TRANSFORMED DATA 
transformed data {                                    // Data block

}


// PARAMETERS PARAMETERS PARAMETERS PARAMETERS PARAMETERS 
parameters {                             // Parameters block
// compound parameters that appear in the probability model see appendix S2
real<lower=0> c2[1]; 
real<lower=0> c3[1];
real<lower=0> c1[1]; 

}

// TRANSFORMS TRANSFORMS TRANSFORMS TRANSFORMS TRANSFORMS
transformed parameters {
  //components of binomial succes
  real pA_T_A[D_NumGrps[1]]; 
  real pA_T_I[1];
  real pB_T_A[1];
  real pB_T_I[D_NumGrps[2]];
  
  //we do separate fitings for AP and IAP because they typically have different time structures
  // AAP AAP AAP AAP AAP AAP AAP AAP AAP AAP AAP AAP AAP
  for (kk in 1:(D_NumGrps[1])) //D_NumGrps
  {
    pA_T_A[kk]=(1-exp(-c3[1]*D_LensA[kk]));
    pB_T_A[1]=(1-exp(-c2[1]*D_bgLens[1,2]));
  }
  // IAP IAP IAP IAP IAP IAP IAP IAP IAP IAP IAP IAP IAP
  for (ss in 1:(D_NumGrps[2])) //D_NumGrps
  {
    pA_T_I[1]=(1-exp(-c3[1]*D_bgLens[2,1]));
    pB_T_I[ss]=(1-exp(-c2[1]*D_LensI[ss]));
  }
  
}

// MODEL MODEL MODEL MODEL MODEL MODEL MODEL MODEL MODEL 
model {     
  
  // half normal prior distributions (becuase compound parameter declarations restrict to non-zero positives)
  c2[1] ~ normal(0,1); 
  c3[1] ~ normal(0,1);  
  c1[1] ~ normal(0,1);  
  
  // Finally the binommial probability model AAP
  for (kk1 in 1:D_NumGrps[1]) //D_NumGrps  
  {
    target += binomial_lpmf(D_InfsA[kk1] | D_RepsA[kk1],1-((1-(c1[1]*pA_T_A[kk1]*pB_T_A[1]))^D_Wf0));
  }
  
  
  // Finally the binommial probability model IAP
  for (uu in 1:D_NumGrps[2])
  {
    target += binomial_lpmf(D_InfsI[uu] | D_RepsI[uu], 1-((1-(c1[1]*pA_T_I[1]*pB_T_I[uu]))^D_Wf0));
  }
  
}
////////////////////////////////////////////////////////////////////  




// GENERATE GENERATE GENERATE GENERATE GENERATE GENERATE 

generated quantities {      // Generated quantities block. 

// extracting alpha, beta and mu from the compuund parameter posteriors by approximation, see appendix S2
real<lower=0> albe[1]; // acquisition
real<lower=0> al[1]; // acquisition
real<lower=0> be[1]; // inoculation
real<lower=0> mu[1]; // virus loss
// and in addition a second round of approximations which may be more accurate but likely ignore this part
real<lower=0> al_p[1]; // acquisition
real<lower=0> be_p[1]; // inoculation
real<lower=0> mu_p[1]; // virus loss

albe[1]=c1[1]*c2[1]*c3[1];
mu[1]=((c3[1]*c2[1])-albe[1])/(c3[1]+c2[1]);
al[1]=c3[1]-mu[1];
be[1]=albe[1]/al[1];

mu_p[1]=((c3[1]*c2[1])-albe[1])/(c3[1]+c2[1]-mu[1]);
al_p[1]=c3[1]-mu_p[1];
be_p[1]=albe[1]/al_p[1];

}

