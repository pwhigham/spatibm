% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixation.R
\name{panmictic.fixation}
\alias{panmictic.fixation}
\title{Time to fixation for a neutral allele}
\usage{
panmictic.fixation(p, N)
}
\arguments{
\item{p}{initial allele frequency}

\item{N}{Population size}
}
\value{
Mean time to fixation
}
\description{
Approximate mean time to fixation for a neutral allele in a well-mixed (panmictic)
population of N individuals with original allele frequency p.

More formally, this is the approximate mean number of generations that a selectively neutral allele initially at a frequency p
will remain segregating in a randomly mating monoecious diploid
population of size N with no spatial structure.  Note that there is no distinction
between males/females and no spatial structure.
This a well-known result in classical population-genetic theory.  See
Ewens, W.J., 1979. Mathematical Population Genetics. Springer-Verlag, Berlin.
For a detailed analysis of fixation due to spatial structure see:
Whigham, Peter A., Dick, Grant C. and Spencer, Hamish G., (2008),
Genetic drift on networks: Ploidy and the time to fixation,
Theoretical Population Biology, 74(4), pp. 283-290,
}
