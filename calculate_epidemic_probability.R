calculate_epidemic_probability <- function(numberInsects,start_intervl,localParameters,virusParameters){
  numVariables <- ((numberInsects+1)*3)-1; 
  
  #### INITIAL section selects the an efficient timestep
  interval_indN=start_intervl
  interval_rec=24*interval_indN
  result <- solveInoculumStatesBP(numVariables, localParameters, virusParameters, (1/interval_rec))
  dr1=dim(result)[1]
  ttilde=result[1:(numVariables+1),]
  f=result[(numVariables+1+1):dr1,]
  # the best choice of interval length is the biggest which produces only positive entries in ttilde - see MS appendix S3
  print('looking for interval length')
  while (sum(ttilde<0)) {
    interval_indN=interval_indN+1
    interval_rec=24*interval_indN
    result=solveInoculumStatesBP(numVariables, localParameters, virusParameters, (1/interval_rec))
    dr1=dim(result)[1]
    ttilde=result[1:(numVariables+1),]
    f=result[(numVariables+1+1):dr1,]
  }
  print('found interval length')
  
  #### SUBSEQUENT section iterates branching process method for finding extinction probability, see MS appendix S3
  itimeCtr <- 0
  k <- dim(ttilde)[2]
  qm <- rep(0, k)           # extinction probability vector
  qprime <- rep(0, k)      # updated extinction probability
  qm_old=rep(-1, k)
  thresh=0.0000001
  # Iterating time loop
  print('start extinction solve')
  while(sum(abs(qm-qm_old)>thresh)>0){#(itimeCtr<100000){#
    qm_old=qm
    itimeCtr=itimeCtr+1
    gb <- 1 - f[numVariables-numberInsects+1, ] + f[numVariables-numberInsects+1, ] * qm[numVariables-numberInsects+1]  # Bernoulli term
    gt <- t(ttilde) %*% c(qm, 1)     # Matrix multiplication (1x5) %*% (5x1)
    qprime <- gt * gb  # update extinction probability
    qm <- qprime  # Update q with qprime for next iteration
  }
  print('end extinction solve')
  
  return(1-qm);
}