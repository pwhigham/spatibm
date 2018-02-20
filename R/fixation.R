#############################################
# Spatial Individual-based Genetic Model
#############################################
# Simple example to model time to fixation of a neutral allele
#
#
# t ( p, N) = -4N  * [ p ln ( p) + (1- p)ln (1- p) ]
#############

#' Time to fixation for a neutral allele
#' @description Approximate mean time to fixation for a neutral allele in a well-mixed (panmictic)
#' population of N individuals with original allele frequency p.
#'
#' More formally, this is the approximate mean number of generations that a selectively neutral allele initially at a frequency p
#' will remain segregating in a randomly mating monoecious diploid
#' population of size N with no spatial structure.  Note that there is no distinction
#' between males/females and no spatial structure.
#' This a well-known result in classical population-genetic theory.  See
#' Ewens, W.J., 1979. Mathematical Population Genetics. Springer-Verlag, Berlin.
#' For a detailed analysis of fixation due to spatial structure see:
#' Whigham, Peter A., Dick, Grant C. and Spencer, Hamish G., (2008),
#' Genetic drift on networks: Ploidy and the time to fixation,
#' Theoretical Population Biology, 74(4), pp. 283-290,
#'
#'@param p initial allele frequency
#'@param N Population size
#'@return  Mean time to fixation
#'
panmictic.fixation <- function(p,N)
{
	-4 * N * ( (p * log(p)) + ((1-p)*log(1-p)))
}

#' Allele fixation
#' @description Has the population lost all allelic diversity?
#' Note that the function uses the last 2 allele values for each individual.
#' Normally fixation would be run for a single loci, so this assumption is
#' not too limiting.
#' @param pop Current population
#' @return TRUE - the population for the final loci has converged or is empty, FALSE otherwise.
allele.fixation <- function(pop)
{
  # Use last two columns (last loci) for alleles

  if (pop$n==0) return(TRUE)

  sum.alleles <- sum(unlist(pop$marks[,(ncol(pop$marks)-1):ncol(pop$marks)]))
	if (sum.alleles==(2*pop$n)) return(TRUE)  # all 1's
	if (sum.alleles==0) return(TRUE)  # all 0's

	FALSE  # mixed
}

#' Time to fixation for an spatial individual-based model
#' @description A simple example of using the \code{single.step} function
#' to determine the time to fixation of a neutral allele for a model.  Assumes
#' that the \code{curr.pop} has been suitably initialised. Note that all parameters for
#' the model, other than \code{max.gens} is the same as the model \code{single.step}.
#' This runs the model once, whereas because of the stochastic nature of the model, you would
#' need to run this model repeatedly to obtain a good mean and variance of the time to fixation.
#' @param max.gens Maximum number of generations before the simulation halts.
#' @param curr.pop The current population as a .ppp
#' @param move.table Movement table as defined by \code{movement.table}.
#' @param survive.table Survival table as defined by \code{survival.table}.
#' @param breed.table Breeding table as defined by \code{breeding.table}
#' @param habitat.surface The habitat (probability) surface. See \code{circle.habitat} for an example, although
#' this surface could be produced in many ways.  The only requirement is that it must cover the entire window
#' associated with the population.
#' @param crowd.table The crowding table as defined by \code{crowding.table}
#' @param crowding.sigma Bandwidth used for density calculation for crowding
#' @param max.dist Maximum distance between breeding parents
#' @param trace.output TRUE - write out generation details to console
#' @return This function will complete when: (a) The alleles have fixed; (b) The maximum number
#' of generations has been reached; or (c) The population has died out. To allow a simple
#' method for distinguishing each of these cases, the model returns a 3 element vector: <TRUE/FALSE fixed>,
#' <generation>, <Final Population Count>.

fixation <- function(max.gens=100,
                             curr.pop,
                             move.table,
                             survive.table,
                             breed.table,
                             habitat.surface=NA,
                             crowd.table=NA,
                             crowding.sigma=0.0,
                             max.dist=5.0,
                             trace.output=FALSE)
{
  if (trace.output)
  {
    alleles <- unlist(curr.pop$marks[,(ncol(curr.pop$marks)-1):ncol(curr.pop$marks)])
    cat(paste("Gen: 0"," N=",curr.pop$n," A:",sum(alleles),"..",sep=""))
    flush.console();
  }
	fixed <- allele.fixation(pop)
	gen <- 1
	while(!fixed)
	{
		curr.pop <- single.step(curr.pop,
		                        move.table,
		                        survive.table,
		                        breed.table,
		                        habitat.surface,
		                        crowd.table,
		                        crowding.sigma,
		                        max.dist)
		if (trace.output)
		{
		  alleles <- unlist(curr.pop$marks[,(ncol(curr.pop$marks)-1):ncol(curr.pop$marks)])
		  cat(paste("Gen:",gen," N=",curr.pop$n," A:",sum(alleles),"..",sep=""))
		  flush.console();
		}
		if (trace.output) if ((gen %% 5)==0) cat("\n")

		fixed <- allele.fixation(curr.pop)
		if ((curr.pop$n==0) | (gen > max.gens) | fixed) break  # Done

		gen <- gen + 1
	}
	alleles <- unlist(curr.pop$marks[,(ncol(curr.pop$marks)-1):ncol(curr.pop$marks)])

	c(fixed,gen,curr.pop$n,length(which(curr.pop$marks$sex==1)),
	  length(which(curr.pop$marks$sex==2)))
}
