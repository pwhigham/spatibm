#############################################
# Spatial Individual-based Model
#############################################
# movement.R
#############
####################################################################
# Operations to move individuals within a timestep
####################################################################

#' Define movement parameters for individuals
#' @description Define the movement parameters for each age class of males and females.  Some simple
#' checking of parameters is done to ensure the correct number of movement parameters are given for the
#' number of age classes.  Note that movement assumes an isotropic space with no preferred direction.
#' @param m.age.class A vector that defines the values for each male age threshold. These are interpreted
#' as less than or equal values. The first value assumes an age from zero to the defined value.
#' @param m.move A list with a two element vector for each age class that defines the movement
#' for this age class.  The movement value is drawn from a normal distribution N(m,sd)
#' @param f.age.class A vector that defines each female age class.
#' @param f.move A list with a two element vector for each age class.
#' @return The movement table. This is passed to the function \code{move}.
#' @examples
#' # Define a small male movement for age class 0 - 5, and large for >5
#' # Female movement is the same for all females.
#' movement.table(m.age.class=c(5,Inf),m.move=list(c(1,0.1),c(5,2)),
#'                f.age.class=Inf,f.move=list(c(4,1)))
movement.table <- function(m.age.class=c(5,10,Inf),
						   m.move=list(c(3,1),c(5,1),c(7,2)),
						   f.age.class=c(3,Inf),
						   f.move=list(c(2,1),c(4,1))
						   )
{
	# Check that age classes correct length for movement params.
	if (length(m.age.class) != length(m.move)) stop("Male age class mismatch with movement params")
	if (length(f.age.class) != length(f.move)) stop("Female age class mismatch with movement params")
	if (length(which(unlist(lapply(m.move,length))!=2))!=0)
	{
		stop("Male movement parameters incorrectly specified - needs 2 per age class")
	}
	if (length(which(unlist(lapply(f.move,length))!=2))!=0)
	{
		stop("Female movement parameters incorrectly specified - needs 2 per age class")
	}
	# Seem ok - although no check that values are valid

	list(m.age.class,m.move,
		 f.age.class,f.move)
}
#' Move the individuals in the population using the defined movement table.
#' @description Moving individuals in the population across space assumes that there is
#' no preferred direction and therefore space is isotropic. Individuals that move outside the
#' defined window for the .ppp population are removed and a warning given.
#' @param p The current population as a .ppp object
#' @param m.table The movement tablee as defined using \code{movement.table}
#' @return A new population with each individual moved according to the \code{movement.table}.
#' Individuals that have moved outside of the window defined by the population are removed and a
#' warning given.  This function should not be directly called by the user.
#
move <- function(p,m.table)
{
	if (p$n == 0) return(p)

	means <- vector(mode="double",length=p$n)
	sds <- vector(mode="double",length=p$n)

	for (i in 1:p$n)
	{
		list.index <- ifelse(p[i]$marks[2]==1,1,3)  # index into table based on sex
		a.c <- m.table[[list.index]] # age classes
		m.p <- m.table[[(list.index+1)]] # movement parameters
		m.index <- which(as.numeric(p[i]$marks[3]) <= a.c)[1] # age
		means[i] <- m.p[[m.index]][1]
		sds[i] <- m.p[[m.index]][2]
	}

	d <- abs(rnorm(p$n,means,sds))

	# d <- abs(rnorm(1,m.p[[m.index]][1],m.p[[m.index]][2]))
	# and now place this randomly in the circle around x at
	# distance d
	rn <- matrix(ncol=2,data=rnorm(2*p$n,mean=0,sd=1))
	mag <- sqrt(apply(rn,1,function(x) sum(x^2)))
	offset <- (rn/mag*d) # offset for each point

	#as.numeric(ind[1:2]) + (rn/mag * d)  # scale to sphere radius

	p$x <- p$x + offset[,1]
	p$y <- p$y + offset[,2]

	#p  - used to just return, but no guarantee that the points lie within the window.
	#
	# Construct a new ppp
	ppp(p$x,p$y,window=p$win,marks=p$marks)  # return the updated points
}

#' Estimate home range for a movement parameter
#' @description Given an initial population and movement table, simulate the movement
#' of the population, starting at the origin, using \code{move.table} for \code{timesteps}.
#' The resulting table gives the convex hull area for each individual in the population.  The purpose
#' of this function is to help relate a estimated organism home range to the movement paramters
#' for the model.  Note that the window for \code{pop} determines the constraints of movement.
#' Individuals that move outside the population window are removed.  The model also commences
#' with the individuals at (0,0), so windows should be centred around the origin.
#'
#' @param  pop The population
#' @param move.table The movement table
#' @param timesteps The number of times movement is applied to each individual.
#' @return Vector of convex-hull areas of length \code{pop$n} representing the movement for each individual.
#'
#'
convex.home.range <- function(pop,move.table,timesteps=5)
{
  pop$x <- rep(0,pop$n)
  pop$y <- rep(0,pop$n) # put at origin

  all.pop <- pop
  for (tstep in 1:timesteps)
  {
    pop <- move(pop,move.table)
    all.pop <- superimpose(all.pop,pop)
  }
  areas <- NULL
  for (i in 1:pop$n)  # for each unique individual determine home range as convex area
  {
    p1 <- all.pop[which(all.pop$marks$id==i),]
    cv <- convexhull(p1)
    areas <- c(areas,area.owin(cv))
  }
  areas
}
