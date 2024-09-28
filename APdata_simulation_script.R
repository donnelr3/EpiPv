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
# THIS script sets up the assay structure and the event rates and calls APdata_simulator which simulates random numbers from exponential distributions


rm(list=setdiff(ls(), "filepath"))

#setwd("C:/Users/donnelr3/OneDrive - University Of Cambridge/Archive_ZG_WaveCassava/rFiles/rStanFiles/RPanalysis")
source("APdata_simulator.R")
source("inoc_durtn_calculator.R")

# set assay structure
nReps=30 # number of reps
numWF=10 # number of wf per rep
# AP exp timeline

# virus rates per hr
al_fix=0.1/60 # acquisition rate
be_fix=1/60 # inoculation rate
gam_fix=(1/2)/60 # latency progression rate
mu_fix=(1/(10))/60 # virus clearance rate

# default durations of acquisition, latent and inoculation periods (to be varied in 1 dimension only in each of the AAP LAP and IAP sub-assays)
T_A=4*60
T_L=2*60
T_I=6*60


# varying in 1 dimension only for each of the AAP LAP and IAP sub-assays
numSlots=11
### AAP ASSAY
T_A_vec=round(seq(0,T_A,length.out=numSlots)[2:numSlots])
inocNum_bywf=matrix(0,numWF,nReps)
TPI_assay1=matrix(0,length(T_A_vec),nReps)
for (nv in 1:length(T_A_vec)) {
  for (nwf in 1:numWF) {
    inocD_byReps=APdata_simulator(c(T_A_vec[nv], T_L, T_I),c(al_fix,be_fix,gam_fix,mu_fix),numWF,1)
    inocNum_bywf[nwf,]=rpois(nReps,lambda=be_fix*inocD_byReps) # rewrite perhaps my own function so it can take a vector of intensities (based on different durations) and return vector of data
  }
  TPI_assay1[nv,]=(colSums(inocNum_bywf)>0)*1
}

#pdf("output_plot1.pdf")
#plot(T_A_vec,rowMeans(TPI_assay1), xlab = "AAP duration, mins", ylab = "Prop. test plants infected",ylim=c(0,1))
#dev.off()  

### LAP ASSAY
T_L_vec=seq(0,T_L,length.out=numSlots)[2:numSlots]
inocNum_bywf=matrix(0,numWF,nReps)
TPI_assay2=matrix(0,length(T_L_vec),nReps)
for (nv in 1:length(T_L_vec)) {
  for (nwf in 1:numWF) {
    inocD_byReps=APdata_simulator(c(T_A, T_L_vec[nv], T_I),c(al_fix,be_fix,gam_fix,mu_fix),numWF,1)
    inocNum_bywf[nwf,]=rpois(nReps,lambda=be_fix*inocD_byReps) # rewrite perhaps my own function so it can take a vector of intensities (based on different durations) and return vector of data
  }
  TPI_assay2[nv,]=(colSums(inocNum_bywf)>0)*1
}

#pdf("output_plot2.pdf")
#plot(T_L_vec,rowMeans(TPI_assay2), xlab = "LAP duration, mins", ylab = "Prop. test plants infected",ylim=c(0,1))
#dev.off()  

### IAP ASSAY
T_I_vec=round(seq(0,T_I,length.out=numSlots)[2:numSlots])
inocNum_bywf=matrix(0,numWF,nReps)
TPI_assay3=matrix(0,length(T_I_vec),nReps)
for (nv in 1:length(T_I_vec)) {
  for (nwf in 1:numWF) {
    inocD_byReps=APdata_simulator(c(T_A, T_L, T_I_vec[nv]),c(al_fix,be_fix,gam_fix,mu_fix),numWF,1)
    inocNum_bywf[nwf,]=rpois(nReps,lambda=be_fix*inocD_byReps) # rewrite perhaps my own function so it can take a vector of intensities (based on different durations) and return vector of data
  }
  TPI_assay3[nv,]=(colSums(inocNum_bywf)>0)*1
}

#pdf("output_plot3.pdf")
#plot(T_I_vec,rowMeans(TPI_assay3), xlab = "IAP duration, mins", ylab = "Prop. test plants infected",ylim=c(0,1))
#dev.off()  


# EXPORT AP assay parameters for later analysis
dput(list(T_A_vec,T_L_vec,T_I_vec,T_A,T_L,T_I,numWF), "assay123_dataList_PT.txt")
# EXPORT AP sub-assay data for later analysis
write.table(TPI_assay1, 'TPI_assay1_PT.dat', append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
write.table(TPI_assay2, 'TPI_assay2_PT.dat', append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
write.table(TPI_assay3, 'TPI_assay3_PT.dat', append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)