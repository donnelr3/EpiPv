\name{estimate_virus_parameters_PT}
\alias{estimate_virus_parameters_PT}
\title{Estimate Virus Parameters for Persistently Transmitted Plant Viruses}
\description{
  This function estimates virus transmission parameters from access period data for a given vector-PT virus-plant combination. 
}
\usage{
  estimate_virus_parameters_PT(assay1, assay2, assay3, offdiag_array, IVD, lsEst, virus)
}
\arguments{
  \item{assay1}{A matrix or data frame containing access period data for the first assay (acquisition access period), structured with the first row as feeding durations, the second row as the number of insect vectors used, and the third row as the number of infections.}
  \item{assay2}{A matrix or data frame containing access period data for the second assay (latent access period), structured similarly to \code{assay1}.}
  \item{assay3}{A matrix or data frame containing access period data for the third assay (inoculation access period), structured similarly to \code{assay1}.}
  \item{offdiag_array}{A numerical array containing the fixed access period durations for those durations that are not varied in the assay.}
  \item{IVD}{An integer value representing the number of insects in a cohort for each replication.}
  \item{lsEst}{An integer value representing the estimated survival of insects in the laboratory (in minutes).}
  \item{virus}{A character string representing the name of the virus for output file naming.}
}
\value{
  Returns a matrix containing the posterior parameter distributions for transmission rates \code{mu}, \code{alpha}, \code{beta}, and \code{gamma}.
}
\details{
  The function runs a Stan model for parameter estimation based on access period assays for persistent transmission viruses. It exports summary parameter estimates and the full Markov chains to the main home directory.
}
\examples{
  # Example of how to call the function
   result <- estimate_virus_parameters_PT(assay1, assay2, assay3, offdiag_array, IVD, lsEst, virus_name)
}
