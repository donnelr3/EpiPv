\name{estimate_virus_parameters_SPT}
\alias{estimate_virus_parameters_SPT}
\title{Estimate Virus Transmission Parameters for SPT}
\description{
    This function estimates virus transmission parameters from access period data for a given vector-SPT virus-plant combination. 
}
\usage{
  estimate_virus_parameters_SPT(assay1, assay2, offdiag_array, IVD, lsEst, virus)
}
\arguments{
  \item{assay1}{A matrix or data frame containing access period data for the first assay (acquisition access period), structured with the first row as feeding durations, the second row as the number of insect vectors used, and the third row as the number of infections.}
  \item{assay2}{A matrix or data frame containing access period data for the second assay (inoculation access period), structured similarly to \code{assay1}.}
  \item{offdiag_array}{A numerical array containing the fixed access period durations for the duration that is not varied in the assay (e.g. where acquisition access period is varied this is typically followed by a fixed inoculationa access period and vice versa).}
  \item{IVD}{An integer value representing the number of insects in a cohort for each replication.}
  \item{lsEst}{An integer value representing the estimated survival of insects in the laboratory (in minutes).}
  \item{virus}{A character string representing the name of the virus for output file naming.}
}
\value{
  Returns a matrix containing the posterior parameter distributions for transmission rates \code{mu}, \code{alpha} and \code{beta}.
}
\details{
  The function runs a Stan model for parameter estimation based on access period assays for persistent transmission viruses. It exports summary parameter estimates and the full Markov chains to the main home directory.
}
\examples{
  # Example usage
   result <- estimate_virus_parameters_SPT(assay1_data, assay2_data, offdiag_data, IVD_value, lsEst_value, "virus_name")
}
