# AP data simulator
# simulates the AAP LAP and IAP assay components of the AP experiment
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
# THIS function simulates the random times of event occurrence

APdata_simulator <- function(lmark_in,smarkpams_in,WF_in,nReps_in) {
  
  WF=WF_in
  nReps=nReps_in
  alrate=smarkpams_in[1]
  berate=smarkpams_in[2]
  gamrate=smarkpams_in[3]
  murate=smarkpams_in[4]
  lmark=lmark_in
  print(lmark[1])
  print(lmark[2])
  print(lmark[3])
  TPI_WFbyRep=matrix(-99,WF,nReps)
  out_durtn=matrix(-99,WF,nReps)
  
  
  for (wf in 1:WF) {
    
    # SIMULATE TIMES UNTIL ACQUISITION 
    tA_sim=rexp(nReps,alrate) # if tA_sim > T_A then no acquisition occurred etc.
    if ((sum(is.na(tA_sim))+sum(is.nan(tA_sim)))>0){
      print('error with the acquisition sim data')
      print(tA_sim)
      print(alrate)
      return(NULL)
    }else{
      #print('confirmed the acquisition sim data as clear of NAs or NANs')
    }
    # SIMULATE TIMES UNTIL LATENCY PROGRESSION 
    tL_sim=rexp(nReps,gamrate)
    if ((sum(is.na(tL_sim))+sum(is.nan(tL_sim)))>0){
      print('error with the latent sim data')
      print(tL_sim)
      print(gamrate)
      return(NULL)
    }else{
      #print('confirmed the latent sim data as clear of NAs or NANs')
    }
    # SIMULATE TIMES UNTIL RECOVERY 
    tR_sim=rexp(nReps,murate)
    if ((sum(is.na(tR_sim))+sum(is.nan(tR_sim)))>0){
      print('error with the retention sim data')
      print(tR_sim)
      print(murate)
      return(NULL)
    }else{
      #print('confirmed the retention sim data as clear of NAs or NANs')
    }

    # provide new acquisition, latency progression, and recovery times if there is recovery before the end of the AAP
    clearEarly=tA_sim+tL_sim+tR_sim
    anyEarly=(clearEarly<lmark[1])
    numReacq=0
    numReacq=numReacq+sum(anyEarly)
    while(sum(anyEarly)>0){
      #resample sum(anyEarly) rv.s
      print(paste('There are ', sum(anyEarly), ' potential re-acquisitions'))
      tA_sim[anyEarly]=rexp(sum(anyEarly),alrate)+clearEarly[anyEarly]
      tL_sim[anyEarly]=rexp(sum(anyEarly),gamrate)
      tR_sim[anyEarly]=rexp(sum(anyEarly),murate)
      clearEarly=tA_sim+tL_sim+tR_sim
      anyEarly=(clearEarly<lmark[1])
      numReacq=numReacq+sum(anyEarly)
    }


    sim_eventTs=cbind(tA_sim,tL_sim,tR_sim)

    tmp_casesDurations=inoc_durtn_calculator(lmark,sim_eventTs) # rewrite so it can work with vectors as there will be simulated nRep based data to pass to this function
    
    if ((sum((tmp_casesDurations<0))>0)||(sum(is.nan(tmp_casesDurations)>0))||(sum(is.na(tmp_casesDurations)>0))){
      print(paste('tmp_casesDurations still contains etiher dummy values of NANs ... tmp_casesDurations dump:'))
      print(tmp_casesDurations)
      return(NULL)
    }
    
    out_durtn[wf,]=tmp_casesDurations[,5]  # IN REALLITY THE SELECTION FROM 1 OF 4 COLUMNS WILL BE A CALCULAITON OR CRITERION
    
    TPI_WFbyRep[wf,]=rpois(nReps,lambda=berate*out_durtn[wf,]) # rewrite perhaps my own function so it can take a vector of intensities (based on different durations) and return vector of data
    # RETURN AND FIX UP NA's where they are caused by 0 inoc durtn
    TPI_WFbyRep[out_durtn[,]==0]=0
    print(paste('There were ', numReacq, ' re-acquisitions for this insect across reps....'))
  }
  
  if ((sum((TPI_WFbyRep<0))>0)||(sum(is.nan(TPI_WFbyRep)>0))||(sum(is.na(TPI_WFbyRep)>0))){
    print(paste('TPI_WFbyRep still contains etiher dummy values of NANs ... TPI_WFbyRep dump:'))
    print(TPI_WFbyRep)
    return(NULL)
  }
  
  any_TPI_WFbyRep=(TPI_WFbyRep>0)*1
  
  TPI=numeric(0)
  
  for (rp in 1:nReps) {
    TPI=c(TPI,(sum(any_TPI_WFbyRep[,rp])>0))
  }

  return(TPI)
}