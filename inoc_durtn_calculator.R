# AP data simulator simulates the AAP LAP and IAP assay components of the AP experiment
# in these experiments we focus on how events occur in these experiments. Namely:
# though it may appear that there is a period for acquiring then a period 
# for passing latency and then a period for inoculating,
# actually once acquisition occurs in the initial period passing 
# latency may occur in the remaining initial period for acquisition, in the 
# specific period for  passing latency, or in the period given for inoculation. 
# When this actually  occurs determines the actual period in which inoculation 
# may occur. As set out in accompanying documentation/text the above 
# consideration leads to a piece-wise inoculation duration depending on when 
# latency progression occurs. The below is a representative simulation of these 
# processes producing an AP dataset given this consideration (actual v apparent)
# THIS function is for calculating the effective period which an insect has for inoculation 

inoc_durtn_calculator <- function(lmark_in,sim_Ts_in) {
  
  ts=dim(sim_Ts_in)
  
  print(paste('size of sim_Ts_in',ts))
  
  tA_sim=sim_Ts_in[,1]
  tL_sim=sim_Ts_in[,2]
  tR_sim=sim_Ts_in[,3]
  
  T_Ain=lmark_in[1]
  T_Lin=lmark_in[2]+T_Ain
  T_Iin=lmark_in[3]+T_Lin
  
  
  output_cases=matrix(-99,ts[1],5) #4 cases - see model - 1-4 are separate cases and 5 is piece-wise. We will only use 5 but 1-4 is useful for checking above against
  
  
  if ((sum(is.na(tA_sim)))>0){
    print('error with little t_A hacing NAs')
    print(tA_sim)
    break;
  }
  if ((sum(is.na(tL_sim)))>0){
    print('error with little t_L having NAs')
    print(tA_sim)
    break;
  }
  if ((sum(is.na(tR_sim)))>0){
    print('error with little t_R having NAs')
    print(tA_sim)
    break;
  }
  
  for (tt in 1:ts[1]) {
    
    
    if (T_Ain<tA_sim[tt]){print('warning no acquisition in acquisition period')}
    if (T_Lin<(tA_sim[tt]+tL_sim[tt])){print('warning no latency in latent period but can still pass in inoc period')}
    if (T_Lin<(tA_sim[tt]+tL_sim[tt]+tR_sim[tt])){print('warning loss of pathogen prior to inoculation period')} 
    if ((tA_sim[tt]>T_Ain)||((tA_sim[tt]+tL_sim[tt])>T_Iin)||((tA_sim[tt]+tL_sim[tt]+tR_sim[tt])<T_Lin)) {
      output_cases[tt,]=rep(0,5)    
      # ALSO want the answer to be no duration if no acquisition occurred BEFORE end of AAP:
      # OR if latency passed AFTER end of IAP
      # OR if retention lost prior to IAP
    }else{
      output_cases[tt,1]=tR_sim[tt]-(T_Lin-(tA_sim[tt]+tL_sim[tt]))
      output_cases[tt,2]=T_Iin-T_Lin
      output_cases[tt,3]=tR_sim[tt]
      output_cases[tt,4]=(T_Iin-T_Lin)-((tA_sim[tt]+tL_sim[tt])-T_Lin)
      output_cases[tt,5]=max((T_Iin-T_Lin)-max(0,(T_Iin-(tA_sim[tt]+tL_sim[tt]+tR_sim[tt])))-max(0,(tA_sim[tt]+tL_sim[tt]-T_Lin)),0)
    }
    
    
    
  }
  
  return(output_cases)
}

