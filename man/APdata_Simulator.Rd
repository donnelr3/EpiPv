\name{APdata_simulator}
\alias{APdata_simulator}
\title{Simulate Access Period Data for Virus Transmission Assays}
\description{
  The \code{APdata_simulator} function receives event rate parameters for virus transmission (\mu, \alpha, \beta, \gamma) as well as assay feeding durations and returns simulated access period data. The simulation considers how events occur during assays, accounting for periods of acquisition, latency progression, inoculation and loss of infectiousness.
}
\usage{
  APdata_simulator(lmark_in, smarkpams_in, WF_in, nReps_in)
}
\arguments{
  \item{lmark_in}{A numeric vector representing assay feeding durations.}
  \item{smarkpams_in}{A numeric vector containing event rate parameters for virus transmission, including \code{alpha}, \code{beta}, \code{gamma} and \code{mu}.}
  \item{WF_in}{An integer indicating the number of whitefly in a cohort per replication}
  \item{nReps_in}{An integer representing the number of replications for the simulation.}
}
\value{
  A binary vector of dimension nReps_in representing the simulated access period data (test plant infection or no test plant infection) based on the input parameters.
}
\details{
  This function simulates the random times of event occurrence in AP (Access Period: Acquisition Access Period, Latent Access Period, and Inoculation Access Period) experiments. The simulation reflects the contingent nature of inoculation duration, which varies depending on when latency progression and loss of infectiousness occurs.
}
\examples{
  # Example of how to use the APdata_simulator function
  lmark <- c(10, 15, 20)  # Assay feeding durations
  smarkpams <- c(0.5, 2, 1, 0.3)  # Event rate parameters: alpha, beta, mu
  WF <- 5  # Number of workflow iterations
  nReps <- 100  # Number of repetitions
  
  # Simulate access period data
  simulated_data <- APdata_simulator(lmark, smarkpams, WF, nReps)
  print(simulated_data)
}
