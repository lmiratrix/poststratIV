% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pack_data_generation_functions.R
\name{rand.exp}
\alias{rand.exp}
\title{Randomize into tx and control and observe outcomes}
\usage{
rand.exp(df, p = 0.3)
}
\arguments{
\item{df}{Dataframe of potential outcomes.  In particular has
columns Y0, Y1, and S0, S1.}

\item{p}{Proportion of units that are treated.}
}
\description{
Useful for running simulation studies.
}
