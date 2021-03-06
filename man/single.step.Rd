% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{single.step}
\alias{single.step}
\title{Perform a single (generation) step for the population}
\usage{
single.step(curr.pop, move.table, survive.table, breed.table, habitat.surface,
  crowd.table = NA, crowding.sigma = 0, max.dist = 5)
}
\arguments{
\item{curr.pop}{The current population as a .ppp}

\item{move.table}{Movement table as defined by \code{movement.table}.}

\item{survive.table}{Survival table as defined by \code{survival.table}.}

\item{breed.table}{Breeding table as defined by \code{breeding.table}}

\item{habitat.surface}{The habitat (probability) surface. See \code{circle.habitat} for an example, although
this surface could be produced in many ways.  The only requirement is that it must cover the entire window
associated with the population.}

\item{crowd.table}{The crowding table as defined by \code{crowding.table}.  May be NA for no crowding.}

\item{crowding.sigma}{Bandwidth used for density calculation for crowding.  Value of 0.0 indicates no crowding.}

\item{max.dist}{Maximum distance between breeding parents}
}
\value{
The next population as a .ppp
}
\description{
Perform a single step for a population based on their defined breeding, survival,
and crowding parameters and the habitat surface.  This function is only called when the
user requires fine-grained control over changes to the population each generation.
}
