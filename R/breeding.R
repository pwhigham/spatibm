#############################################
# Spatial Individual-based (Genetic) Model
#############################################
# breeding.R
#############
library(spatstat)

#' Define breeding parameters
#' @description Breeding occurs between males and females who are spatially located
#' near each other.  The probability of breeding depends on the age (which may be different
#' for males and females).  Given a breeding event, fecundity is drawn from a defined normal
#' distribution.  Children are randomly placed around the mid-point between parents.
#'
#' @param m.age.class A vector that defines the values for each male age threshold. These are interpreted
#' as less than or equal values. The last value should
#' be Inf.  An interval (lower,upper) is interpreted as lower <= V <= upper.  The
#' first value is assumed to be the upper value for a range (0,upper).
#' @param m.breed A vector of probabilities of breeding, one for each male age class
#' @param f.age.class A vector that defines each female age class.
#' @param f.breed A vector of probabilities of breeding, one for each female age class
#'@param fecund Number of children created if breeding occurs. Defined as two value vector
#'N(m,sd) drawn from a normal distribution.  Values always >= 0.
#'@param displace Displacement of children drawn from normal distribution around the mid-point in space
#'between the male and female parents.
#' @return Breeding table to be used in \code{baseline.breeding}.  Note some basic checking of age class
#' matching to number of probabilities is done, but no checking of values.
#'
###################################################################
# Defining baseline breeding parameters based on age class
###################################################################
breeding.table <- function(m.age.class=c(5,40,Inf),
                           m.breed=c(0.0,0.8,0.1),
                           f.age.class=c(6,40,Inf),
                           f.breed=c(0.8,0.9,0.1),
                           fecund=c(5,2),  # Number of childer as N(m,sd)
                           displace=c(3,0.5) # Displacement from mid-point of parents as N(m,sd)
                           )
{
	# Check that age classes correct length for movement params.
	if (length(m.age.class) != length(m.breed)) stop("Male age class mismatch with breeding params")
	if (length(f.age.class) != length(f.breed)) stop("Female age class mismatch with breeding params")

  	# Seem ok - although we haven't checked values make any sense!

	list(c(0,m.age.class),m.breed,
		 c(0,f.age.class),f.breed,
		 fecund,displace)
}
#' Create population available for breeding
#' @description Given a population and a breeding table, determine which individuals
#' can breed based on their age class and the probability of breeding.
#'@param p Current population of individuals as .ppp
#'@param b.table The breeding table defined by \code{breeding.table}
#'@return A new population containing only those individuals who are available for breeding.
#
baseline.breeding <- function(p,b.table)
{
	if (p$n==0) return(p)
	# For each individual, determine if they can breed
	#
	breed.surv <- vector(mode="double",length=p$n)

	males <- which(p$marks$sex==1)
	females <- which(p$marks$sex==2)

	if (length(males) > 0)
	{
    i <- cut(p$marks$age[males],breaks=b.table[[1]],include.lowest=TRUE)
    breed.surv[males] <- unlist(b.table[[2]][i])
	}
	if (length(females) > 0)
	{
	  i <- cut(p$marks$age[females],breaks=b.table[[3]],include.lowest=TRUE)
	  breed.surv[females] <- unlist(b.table[[4]][i])
	}
	p[(runif(p$n) < breed.surv),]
}

#' Construct breeding pairs of males/females
#' @description Given a set of individuals who are available for breeding, pair
#' up the males/females.  This is done by finding the first female and matching to closest (geographic)
#' male.  This pair is then removed from the available pool for breeding, and the process
#' repeated.  Valid breeding pairs must be no further apart than \code{max.dist}.
#'@param p The population of valid breeding individuals
#'@param max.dist Maximum distance allowed between male/female pairs for breeding
#'@return A 2 column table (female, male) of individual indices of breeding pairs.
#'
breeding.pairs <- function(p,max.dist)
{
	if (p$n == 0) return(NULL)

	b.pairs <- NULL # Breeding pairs to collect as a table
	#
	pair.d <- pairdist(p) # Get all distances between individuals
	                 # Square matrix nrow,ncol= p$n
	# Make matrix with rows as females, columns as males
	f <- which(p$marks$sex==2)
	m <- which(p$marks$sex==1)

	while((length(f) > 0) & (length(m) > 0))
	{
		d <- pair.d[f,m]
	# Now select the female with the closest males out of all females
		min.d <- min(d)

		if (min.d > max.dist) break # we are finished

		fem <- which(apply(matrix(d,ncol=length(m)),1,min)==min.d)
		mal <- which(apply(matrix(d,ncol=length(m)),2,min)==min.d)
		b.pairs <- rbind(b.pairs,c(f[fem],m[mal]))

	# and remove them from the list of f/m candidates for breeding
	#
		f <- f[-fem]
		m <- m[-mal]
	}
	return(b.pairs)
}
#' Breed a pair of individuals
#' @description Breed a pair of individuals, creating a set of offspring.  These are
#' placed in space based on \code{displace}.  Medellian genetics are used to determine the
#' new allele values for the offspring.
#' @param x A row from the breeding pairs table, giving the male and female ids.
#' @param pts The .ppp population of breeding individuals.  This is used to access the (x,y) location
#' and allele values for the parents.
#'@param fecund Number of children created if breeding occurs. Defined as two value vector
#'N(m,sd) drawn from a normal distribution.  Values always >= 0.
#'@param displace Displacement of children drawn from normal distribution around the mid-point in space
#'between the male and female parents.
#' @return The table of children.  See \code{breed} for details.
##########################
do.breeding <- function(x,pts,fecund,displace)
{
	xy <- c(((pts$x[x[1]]+pts$x[x[2]])/2),
			((pts$y[x[1]]+pts$y[x[2]])/2))

	num.children <- floor(rnorm(1,mean=fecund[1],sd=fecund[2]))

	if (num.children <= 0) return(NULL)

	res <- matrix(nrow=num.children,ncol=(ncol(pts$marks)+2))
	colnames(res) <- c("x","y",colnames(pts$marks))



	res[,4] <- ifelse(runif(num.children) > 0.5,1,2)
	rn <- matrix(ncol=2,data=rnorm(2*num.children,0,1)^2)
	mag <- sqrt(rowSums(rn))
	d <- abs(rnorm(num.children,displace[1],displace[2]))
	res[,1:2] <- matrix(ncol=2,data=(xy + (rn/mag * d)),byrow=TRUE)
	res[,6] <- pts$marks[x[1],1]
	res[,7] <- pts$marks[x[2],1]

	#########################################
	# MENDELIAN Inheritance for diploid
	#########################################

	# Make combined table of parent alleles
	palleles <- cbind(as.numeric(pts$marks[x[1],6:ncol(pts$marks)]),
	                  as.numeric(pts$marks[x[2],6:ncol(pts$marks)]))

 # seln <- sample(1:2,(ncol(pts$marks)-5)*num.children,replace=T)

	res[,8:ncol(res)] <- t(replicate(num.children, palleles[cbind(seq(nrow(palleles)), sample(1:2, nrow(palleles), replace=TRUE))]))

	res
}
#' Breed valid individuals
#' @description For each closest pair of male/females, create \code{fecund}
#'  number of children, and place them
#'     around the midpoint of the location between the parents, displaced
#'     from this point by \code{displace} as an isotropic placement.
#'     Once parents are seperated by > \code{max.dist} breeding stops.
#' @param pts The points representing individuals that are available for
#'           breeding
#'@param b.pairs The breeding pairs as a matrix returned from \code{breeding.pairs}.
#'      Each row is <female> <male> index into \code{pts}.
#'@param b.table Breeding table.  This is used to extract the number of children to produce and
#' displacement of children drawn from normal distribution around the mid-point in space
#' between the male and female parents.
#' @return  The children as a .ppp with id,sex,age, allele values updated.
#############################################################################
breed <- function(pts,b.pairs,b.table)
{
  fecund <- b.table[[5]]
  displace <- b.table[[6]]
	kids <- NULL
	for (parents in 1:nrow(b.pairs))
	{
		kids <- rbind(kids,
					  do.breeding(b.pairs[parents,],
							pts,fecund,displace))
	}
	if (is.null(kids)) return(NULL)  # Failed to breed

	# Update all kids age and ids

	kids[,5] <- 0 # age
	kids[,3] <- 0 # set id as 0 - updated in multi-step

	kids.marks <- as.data.frame(kids[,3:ncol(kids)])

	if (ncol(kids.marks)==1) kids.marks <- t(kids.marks)  # hack for R

	ppp(x=kids[,1],y=kids[,2],
		window=pts$window,
		marks=kids.marks)  # Return the point dataset of children
}





