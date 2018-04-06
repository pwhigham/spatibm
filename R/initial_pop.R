#############################################
# Spatial Individual-based Model
#############################################
# initial_pop.R
################################################################
library(spatstat)

# Creating initial population and spatial description
# To generate documentation:  devtools::document()
################################################################

################################################################################
#' Generate a Poisson stationary set of N points within the observation
#' window owin.  The intensity of the points is given by N / area(owin)
#' @param N Number of points in the population
#' @param owin The observation window
#' @return A .ppp point pattern based on the rpoispp (spatstat) function
#' @examples
#' pts <- make_rpois_num(5,
#'                       owin=owin(xrange=c(-10,10),
#'                                 yrange=c(-10,10),unitname="metre"))
################################################################################
make_rpois_num <- function(N,owin)
{
	N.intensity <- N/area.owin(owin)
	while(TRUE)
	{
		rec.sim <- rpoispp(N.intensity,win=owin)
		if (rec.sim$n==N) return(rec.sim)
	}
}

########################################################################################################
# create.ibm.population
#' Create an initial population
#' @description Create an initial population by defining a release window, spatial layout, female/male distribution,
#' age distribution and allele proportions for each loci
#'
#' @param N Number of individuals in population
#' @param spat.layout String to select layout options.  Currently only
#'                     "random" is implemented.
#' @param obs.owin  Observation window that defines the boundaries of the space where
#'                  individuals are placed.
#' @param prob.females  Probability of females in population (0 - 1)
#' @param age.distribution The age distribution for the population.  Drawn from
#'                        a normal distribution c(mean, sd) Age must be >= 1.
#'                        Generated age is an integer. Separate age distribution properties may be
#'                        given for females by setting age.distribution to a four value vector, in which case
#'                        age distribution is c(male.mean, male.sd, female.mean,female.sd)
#' @param allele.prop  Proportion of allele Ai per ith loci for the population.
#'                     Well mixed = 0.5, converged = 1.0 (or 0.0)
#'                     Alleles are stored in the data frame as a vector of size 2 * number loci.
#'                     This is done to make the breeding cycle easier when it is modelled.
#'                     For example, with 2 alleles would be represented as (A1,B1,A2,B2) where A1 and A2
#'                     are the values on the 2 loci first strand, B1 B2 matching loci on second strand.
#'
#' @return A .ppp object containing N points and a marks dataframe that defines
#' for each point a unique id, sex, age, male parent id, female parent id, and allele values.
#'
#' @examples
#' # Create the default population and plot using plot.ppp
#' curr.pop <- create.ibm.population()
#' plot(curr.pop,use.marks=F,main="Initial Pop")

create.ibm.population <- function(N=50,
                                  spat.layout="random",
					                        obs.win=owin(xrange=c(-10,10),
					                                     yrange=c(-10,10),unitname="metre"),
					                        prob.females=0.5,
					                        age.distribution=c(10,1,7,2),
					                        allele.prop=c(0.5,0.9))
{
	# SEX
	f = runif(N)
  sex <- rep(1,N) # Males = 1
	sex[f <= prob.females] <- 2 # Females = 2
	#
	# ID
	id <- seq(from=1,to=N)
	# Parent IDs - set to 0 since all IDS >= 1
	m <- rep(0,N)
	f <- rep(0,N)
	# AGE
	age <- ceiling(abs(rnorm(N,mean=age.distribution[1], sd=age.distribution[2])))
  if (length(age.distribution) > 2) # Separate values for females
  {
    age[which(sex==2)] <- ceiling(abs(rnorm(length(which(sex==2)),mean=age.distribution[3], sd=age.distribution[4])))
  }
	age[age<=0] <- 1  # make sure age at least 1
    # ALLELES
	alleles <- as.data.frame(matrix(nrow=N,ncol=2*length(allele.prop),data=1))
	index <- 1
	for (p in allele.prop)
	{
		alleles[(runif(N) < p),index] <- 0
		alleles[(runif(N) < p),(index+1)] <- 0
		index <- index + 2
	}
	pop_data <- as.data.frame(cbind(id,sex,age,f,m,alleles))  # LOSES sex as a factor ????
	#
	# Now just generate the points and connect the data....
	#
	pop <- switch(spat.layout,
				"random" = make_rpois_num(N,obs.win)
		   )
	setmarks(pop,pop_data)
}
