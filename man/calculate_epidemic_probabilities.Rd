\name{calculate_epidemic_probability}
\alias{calculate_epidemic_probability}
\title{Calculate Epidemic Probability given Plant Virus Transmission rates and local parameters}
\description{
  This function calculates the epidemic probability for different types of inoculum state based on event rate parameters for virus transmission, local parameters, and user inputs.
}
\usage{
  calculate_epidemic_probability(numberInsects, interv, localParameters, virusParameters)
}
\arguments{
  \item{numberInsects}{An integer representing the local level for insect vector burden i.e. the number of vectors per plant.}
  \item{interv}{A numeric value indicating the initial time interval resolution (minutes) for the probability calculations (recommend using 10 and the routine will increase this by 1 if step size is too big until step size is appropriate (i.e., probability no event occurs >0)).}
  \item{localParameters}{A numeric vector containing local event rate parameters (\theta, r, h, b_f, nu_pl) (vector dispersal, roguing, harvesting, vector field mortality, plant latency progression rates, all day^{-1}).}
  \item{virusParameters}{A numeric vector containing virus transmission parameters (\alpha, \beta, \mu) (acquisition, inoculation and virus clearance rates, all day^{-1}).}
}
\value{
  A numeric array or matrix representing the epidemic probabilities for various inoculum states.
}
\examples{
  # Example of how to use the calculate_epidemic_probability function
  numberInsects <- 10
  interv <- 10
  localParameters <- c(1, 0.1, 0.01, 0.1, 0.05)
  virusParameters <- c(0.3, 0.4, 0.5)
  
  epidemic_probabilities <- calculate_epidemic_probability(numVars,starting_guess, numberInsects, interv, localParameters, virusParameters)
  print(epidemic_probabilities)
}
