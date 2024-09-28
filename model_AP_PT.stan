// Description: Probability model of PT access period assay used in Bayesian analysis, see Donnelly Tankam Gilligan 2024 appendix S1 
// data : access period data NOTE that we assume there are acquisition access latent access and inoculation access components for PT virus
//
// infer alpha beta gamma mu from data 
// note that mu in the first instance is rate of loss of insect infectiosness
// in the PT model this includes insect lab mortality
// therefore we additionally use user-inputted expected lab insect survival to separate out virus clearance
// note convention throughout is that variables with D_ in front all involve data coming as input


functions {
  real special_taylor(real pivot, real node, int t_in, int precis_in, int orderExpand){
    real tmp;
    real output;
    if (fabs(node-pivot)<(10^(-precis_in))){
      tmp=0.0;
      for (ss in 1:orderExpand){
        tmp=tmp+((((pivot-node)*t_in)^(ss-1))/tgamma(ss+1));
      }
      output=node*t_in*tmp;
    }else{
      output=(node/(node-pivot))*(1-exp(-(node-pivot)*t_in));
    }
    return output;
  }
  real geom_series(real pivot2_node2, int precis_in, int orderExpand){
    real tmp;
    real output;
    if (fabs(pivot2_node2-1)<(10^(-precis_in))){
      tmp=0.0;
      for (ss in 1:orderExpand){
        tmp=tmp+((pivot2_node2)^(ss-1));
      }
      output=tmp;
    }else{
      output=1/(1-(pivot2_node2));
    }
    return output;
  }
}

// DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA 
data {                                    // Data block
int<lower=0> D_NumGrps[3];  // number of groups to be studied in each of AAP LAP IAP variants (assumes symmetry in exp structure)
int<lower=1> D_Wf0;        // ultimate number of whitefly ON TEST PLANTS
real<lower=1> D_Ls;        // expected insect lifespan in minutes
int<lower=0> D_LensA[D_NumGrps[1]];        // Number of reps for each of AAP LAP IAP variants across studied groups
int<lower=0> D_LensL[D_NumGrps[2]];
int<lower=0> D_LensI[D_NumGrps[3]];
int<lower=0> D_RepsA[D_NumGrps[1]];        // Number of reps for each of AAP LAP IAP variants across studied groups
int<lower=0> D_RepsL[D_NumGrps[2]]; 
int<lower=0> D_RepsI[D_NumGrps[3]]; 
int<lower=0> D_InfsA[D_NumGrps[1]];                        // Response variable (if transm)
int<lower=0> D_InfsL[D_NumGrps[2]];  
int<lower=0> D_InfsI[D_NumGrps[3]];  
int<lower=0> D_bgLens[3,3];
}

// TRANSFORMED DATA TRANSFORMED DATA TRANSFORMED DATA 
transformed data {                                    // Data block

}


// PARAMETERS PARAMETERS PARAMETERS PARAMETERS PARAMETERS 
parameters {                             // Parameters block
real<lower=0,upper=1> mu[1]; 
real<lower=0> al[1]; 
// several parameters are compound for smoother estimation
real<lower=0> mu_lat[1];
real<lower=0> mu_be[1];
}


// TRANSFORMS TRANSFORMS TRANSFORMS TRANSFORMS TRANSFORMS
transformed parameters {
  real<lower=0> lat[1]; // latent progress rate 
  real<lower=0> be[1]; // inoculation rate
  real<lower=0,upper=1> binPam_totsA[D_NumGrps[1]];
  real<lower=0,upper=1> binPam_totsL[D_NumGrps[2]];
  real<lower=0,upper=1> binPam_totsI[D_NumGrps[3]];
  real<lower=0,upper=1> binPam_succA[D_NumGrps[1]];
  real<lower=0,upper=1> binPam_succL[D_NumGrps[2]];
  real<lower=0,upper=1> binPam_succI[D_NumGrps[3]];
  
  lat[1]=mu[1]/mu_lat[1];
  be[1]=mu[1]/mu_be[1];
  
  binPam_totsA = rep_array(0.0,D_NumGrps[1]); // note that stage 0 is handled separately in model fitting code 
  binPam_totsL = rep_array(0.0,D_NumGrps[2]); // note that stage 0 is handled separately in model fitting code 
  binPam_totsI = rep_array(0.0,D_NumGrps[3]); // note that stage 0 is handled separately in model fitting code 
  binPam_succA = rep_array(0.0,D_NumGrps[1]); // note that stage 0 is handled separately in model fitting code 
  binPam_succL = rep_array(0.0,D_NumGrps[2]); // note that stage 0 is handled separately in model fitting code 
  binPam_succI = rep_array(0.0,D_NumGrps[3]); // note that stage 0 is handled separately in model fitting code 
  
  {
    real termA1;
    real termA2;
    real termA3;
    real termB;
    real temp1;
    real temp2;
    real temp3;
    real temp4;
    int TAAP_assay;
    int TLAP_assay;
    int TIAP_assay;
    int orderOfExpand=6;
    int precisForSing=10;
    int TAAP_assay_vec[3];
    int TLAP_assay_vec[3];
    int TIAP_assay_vec[3];
    
    for (dd in 1:(3)) //assays
    {
      TAAP_assay=TAAP_assay_vec[dd];
      TLAP_assay=TLAP_assay_vec[dd];
      TIAP_assay=TIAP_assay_vec[dd];    
      
      for (kk in 1:(D_NumGrps[dd])) //D_NumGrps
      {
        // RECHECK!
        if (dd==1){
          TAAP_assay = D_LensA[kk];
          TLAP_assay = D_LensA[kk]+D_bgLens[1,2];
          TIAP_assay = D_LensA[kk]+D_bgLens[1,2]+D_bgLens[1,3];
        }else if (dd==2){
          TAAP_assay = D_bgLens[2,1];
          TLAP_assay = D_bgLens[2,1]+D_LensL[kk];
          TIAP_assay = D_bgLens[2,1]+D_LensL[kk]+D_bgLens[2,3];
        }else{
          TAAP_assay = D_bgLens[3,1];
          TLAP_assay = D_bgLens[3,1]+D_bgLens[3,2];
          TIAP_assay = D_bgLens[3,1]+D_bgLens[3,2]+D_LensI[kk];
        }
        
            temp1=special_taylor(lat[1],al[1], TAAP_assay, precisForSing, orderOfExpand);
            temp2=special_taylor(mu[1],al[1], TAAP_assay, precisForSing, orderOfExpand);
            temp3=special_taylor((mu[1]+be[1]),lat[1], (TLAP_assay-TIAP_assay), precisForSing, orderOfExpand);

            temp4=geom_series(mu_lat[1], 30, orderOfExpand);

            termA1=temp1*(exp(-(lat[1])*(TLAP_assay)));
            
            termA2=temp2*(exp(-(mu[1])*(TLAP_assay)));
            
            termA3=temp4*(1/(mu_be[1]+1))*(exp(-(mu[1]+be[1])*(TIAP_assay-TLAP_assay))-1);
            
            termB=temp1*(1/(mu_be[1]+1))*(exp(-(lat[1])*(TIAP_assay)))*(temp3-(1-exp(-(lat[1])*(TLAP_assay-TIAP_assay))));
            
            if (dd==1){
              binPam_totsA[kk]=((termA1-termA2)*termA3)+termB;
              binPam_succA[kk]=1-((1-((binPam_totsA[kk])))^D_Wf0);
            }else if (dd==2){
              binPam_totsL[kk]=((termA1-termA2)*termA3)+termB;
              binPam_succL[kk]=1-((1-((binPam_totsL[kk])))^D_Wf0);
            }else{
              binPam_totsI[kk]=((termA1-termA2)*termA3)+termB;
              binPam_succI[kk]=1-((1-((binPam_totsI[kk])))^D_Wf0);
            }
        }
    }
  }
  
}

// MODEL MODEL MODEL MODEL MODEL MODEL MODEL MODEL MODEL 
model {     
  
  
  // cauchy distributed priors
  mu[1] ~ cauchy(1, 1);   
  mu_lat[1] ~ cauchy(1, 1);
  al[1] ~ cauchy(1, 1);   
  mu_be[1] ~ cauchy(1, 1);
  
  for (kk in 1:D_NumGrps[1]) //D_NumGrps
  {
    // AAP variants
    target += binomial_lpmf(D_InfsA[kk] | D_RepsA[kk], binPam_succA[kk]);
  }
  
  for (kk in 1:D_NumGrps[2]) //D_NumGrps
  {
    // LAP variants
    target += binomial_lpmf(D_InfsL[kk] | D_RepsL[kk], binPam_succL[kk]);
  }
  
  for (kk in 1:D_NumGrps[3])
  {
    // IAP variants
    target += binomial_lpmf(D_InfsI[kk] | D_RepsI[kk], binPam_succI[kk]);
  }
  
}
////////////////////////////////////////////////////////////////////  




// GENERATE GENERATE GENERATE GENERATE GENERATE GENERATE 

generated quantities {      // Generated quantities block. 

real mu_cl[1]; // separated out virus clearance rate

mu_cl[1]=mu[1]-(1/D_Ls);


}

