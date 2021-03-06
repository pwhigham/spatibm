% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualise.density.R
\name{origin.density}
\alias{origin.density}
\title{Create a table of densities about the origin}
\usage{
origin.density(win, num.pts = 20, trials = 10, sigma = 1)
}
\arguments{
\item{win}{Spatial window of type \code{owin}}

\item{num.pts}{Maximum number of points to randomly place in the window}

\item{trials}{Number of times to resample points}

\item{sigma}{The bandwidth used for \code{density.ppp}}
}
\value{
Matrix with number of rows as points increase from \code{1:num.pts}, each column
is the measured density at (0,0) for the particular number of points for each trial.
}
\description{
For a given window, randomly place \code{1:num.pts} in a window
and measure the density at (0,0) using the bandwidth \code{sigma}.  This can be
used to understand the density for randomly distributed populations when
setting the crowding parameters for the model.
}
\examples{
# Create a window centred on (0,0)
win <- owin(xrange=c(-2,2),yrange=c(-2,2),unitname="metres")
den <- origin.density(win,num.pts=20,trials=50)
boxplot(t(den),xlab="Number Points",ylab="Density")
}
