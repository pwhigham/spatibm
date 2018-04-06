#############################################
# Spatial Individual-based (Genetic) Model
#############################################
# survival.R
#############
library(spatstat)
library(raster)
####################################################################
# Example operations to handle survival to next breeding
####################################################################

#' Define survival parameters for individuals
#' @description Define the survival parameters for each age class of males and females.  Some simple
#' checking of parameters is done to ensure the correct number of parameters are given for the
#' number of age classes.  Note that survival is measured as a probability (although possibly modified by the habitat).
#' @param m.age.class A vector that defines each male age class.
#' An interval (lower,upper) is interpreted as lower <= V <= upper. The first value is assumed
#' to have a lower bound of zero.
#' @param m.survive A vector of probabilities of survival,  one for each age class
#' @param f.age.class A vector that defines each female age class.
#' @param f.survive A vector of probabilities, one for each female age class.
#' @return The survival table. This is passed to the function \code{baseline.survival}.
survival.table <- function(m.age.class=c(5,40,Inf),
						   m.survive=c(0.8,0.8,0.1),
						   f.age.class=c(3,40,Inf),
						   f.survive=c(0.8,0.9,0.1)
						   )
{
	# Check that age classes correct length for movement params.
	if (length(m.age.class) != length(m.survive)) stop("Male age class mismatch with survival params")
	if (length(f.age.class) != length(f.survive)) stop("Female age class mismatch with survival params")

  # Seem ok - although we haven't checked values make any sense!

	list(c(0,m.age.class),m.survive,
		 c(0,f.age.class),f.survive)
}

#' Baseline survival for individuals
#' @description Given a population \code{p}, determine the survival to the next generation based on
#' the baseline probabilities of survival, \code{s.table}, and the habitat.
#' @param p The current population
#' @param s.table The survival probability table
#' @param habitat The habitat surface representing the probability of survival at any location.
#' @return The surviving population as a .ppp
#'

baseline.survival <- function(p,s.table,habitat=NA)
{
	if (p$n==0) return(p)

	age.surv <- vector(mode="double",length=p$n)

	males <- which(p$marks$sex==1)
	females <- which(p$marks$sex==2)

	if (length(males) > 0)
	{
	  i <- cut(p$marks$age[males],breaks=s.table[[1]],include.lowest=TRUE)
	  age.surv[males] <- s.table[[2]][i]
	}
	if (length(females) > 0)
	{
	  i <- cut(p$marks$age[females],breaks=s.table[[3]],include.lowest=TRUE)
	  age.surv[females] <- s.table[[4]][i]
	}

		# and now do the habitat surface
	if (length(habitat) < 2)
	{
	  final.prob <- age.surv
	}
	else
	{
		hab.prob <- raster(habitat)[cellFromXY(raster(habitat),cbind(p$x,p$y))] # Adjust x,y for raster model
		final.prob <- age.surv*hab.prob
	}

	p[(runif(p$n) < final.prob),]
}

# Plot one of male/female
# ... other plotting parameters
plot.ct <- function(m.t,main="Males",...)
{
	big <- which(m.t$age==Inf)
	max.age <- max(m.t[-big,]$age)  # if only Inf then this is -Inf
	if (max.age==-Inf) max.age <- 0
	m.t[big,"age"] <- max.age + 1
	big.d <- which(m.t$density==Inf)
	max.den <- max(m.t[-big.d,]$density) # ditto
	if (max.den==-Inf) max.den <- 0
	m.t[big.d,"density"] <- max.den + 1

	# Make inf age coloured red...
	cols <- rep(1,nrow(m.t))
	cols[big] <- 2
	syms <- rep(19,nrow(m.t))

	plot(x=m.t$age,y=m.t$density,pch=syms, col=cols,
		cex=m.t$prob*3,xlab="Age",
		ylab="Crowding Density",
		ylim=c(0,max(m.t$density)+2),xlim=c(0.5,max(m.t$age)),main=main,...)
	text(x=m.t$age-0.3,y=m.t$density,labels=m.t$prob,cex=0.7)
	if (length(big.d) > 0)
	{
		abline(h=(max.den+1),lty=2,lwd=0.5,col='grey');
		text(x=1,y=(max.den+1.5),"Infinite Density",cex=0.6)
	}
	if (length(big) > 0) abline(v=(max.age+1),lty=2,lwd=0.5,col='grey')
}
#########################################################################
#' Plot the crowding table
#' @description Graphically show the crowding table as a plot of age versus crowding density.  The probability
#' of survival based on density is shown by value and size of symbol.
#' @param c.t Defined crowding table from \code{crowding.table}.
#' @param ... Graphical parameters
#' @examples
#'ct <- crowding.table()  # Use defaults
#'plot.crowding.table(ct)

plot.crowding.table <- function(c.t,sex=1,...)
{
	plot.ct(subset(c.t,sex==sex),...)
}

make.ctable.entries <- function(sex=1,age.class,den.class,prob)
{
	c.t <- as.data.frame(matrix(ncol=4,nrow=length(prob)))
	colnames(c.t) <- c("sex","age","density","prob")

	c.t$sex <- sex
	i <- 1  # Index into prob
	for (age in age.class)
	{
		for (den in den.class)
		{
			c.t$age[i] <- age
			c.t$density[i] <- den
			c.t$prob[i] <- prob[i]
			i <- i + 1 # Next record
		}
	}
	c.t
}

#################################################################################################
# crowding.table
#################

# IN: m.age.class: age classes for males
#     m.den.class: density classes for each age class
#     m.prob:  probabilities for each age/density combination.  Taken in order of age class.
#     f.age.class, etc.  Female definition.
#
# OP: Checks that basic definition matches number of probabilities per den/age class.
#     and builds the final crowding table.
#
# OUT: The crowding table
#
# The resulting table is passed to density.survival and can be used to visualise the
# crowding definition using plot.crowding.table(c.table)
#
###############################################################################################
#' Define crowding table
#' @description Define the crowding table for the population
#' This involves defining (for both males and females) the age classes and then
#' for each age class, the density classes, and for each age class/density Class
#' the probability of survival given this density measure for the specific age class.
#' To make the definition somewhat simplier we assume that the number of density classes
#' are the same for each age class.  This isn't a limitation, and makes it somewhat easier to
#' define the probabilities, since you just have to have the same number per age class.
#' To make it simplier (and conceptually easier to think about) we split the male and female definitions.
#'@param m.age.class Male age classes
#'@param m.den.class Male density classes
#'@param m.prob Male probabilities of survival given any age and density class
#'@param f.age.class Female age classes
#'@param fden.class Female density classes
#'@param f.prob Female probabilities of survival given any age and density class
#'@return The crowding table.  Note some basic error checking is performed.  The resulting table should be visualised
#'using \code{plot.ct} to confirm the crowding table is correctly defined.
#'@examples
#'ct <- crowding.table()  # Use defaults
#'plot.crowding.table(ct)
#'
crowding.table <- function(m.age.class=c(5,Inf), # age classes
						   m.den.class=c(0.5,1,4,Inf), # same number of density classes for each age
						   m.prob=c(1,0.9,0.5,0.01, # probabilities survival per density class per age
								    1,0.9,0.5,0.01),
						   f.age.class=c(4,Inf), # age classes
						   f.den.class=c(0.5,1,2,4,Inf), # same number of density classes for each age
						   f.prob=c(1,0.9,0.8,0.1,0.01, # probabilities survival per density class per age
								    1,0.9,0.8,0.1,0.01))
{
# Check definition fits constraints.
	if ((length(m.age.class)*length(m.den.class)) != length(m.prob))
	{
		stop("Crowding table: Male probabilities do not match age * den classes")
	}
	if ((length(f.age.class)*length(f.den.class)) != length(f.prob))
	{
		stop("Crowding table: Female probabilities do not match age * den classes")
	}
	c.t <- make.ctable.entries(sex=1,m.age.class,m.den.class,m.prob)
	rbind(c.t,make.ctable.entries(sex=2,f.age.class,f.den.class,f.prob))
}
#' Determine survival due to crowding
#' @description For each individual, determine the density of the remaining population
#' around this individual and then determine the probability of survival as defined in \code{c.table}.
#' Use this to calculate survival. Density is determined by considering the density around an individual, but
#' not including this individual in the density calculation.  The bandwidth of the kernel used for the density
#' calculation is given by \code{sigma}.
#' @param p The current surviving population
#' @param c.table The crowding (density) table, as defined by \code{crowding.table}.
#' @param sigma  The kernel density bandwidth
#' @return The surviving population as a .ppp

crowding.survival <- function(p,c.table,sigma=1)
{
# leaveoneout=TRUE means you don't include the point around which you
# are estimating density...this makes sense since we want to assess the
# density of the neighbourhood.
# sigma is the standard deviation of the isotropic Gaussian kernel
#
	if (p$n==0) return(p)

	d.map <- density.ppp(p,sigma=sigma,at="points",leaveoneout=TRUE)
#
# Ok - for each individual, calculate the probability of survival based on their measured density
#
	res <- rep(1,p$n)
	for (i in 1:p$n)
	{
		subc <- subset(c.table,(c.table$sex==p$marks$sex[i]) &
						(p$marks$age[i] <= c.table$age) &
						(d.map[i] <= c.table$density))
		# First entry is correct one to use
		res[i] <- subc[1,"prob"]
	}
	p[(runif(p$n) <= res),]
}

#' Create plot showing probability of survival over time
#' @description A simple plot showing how the probability of survival reduces
#' over the years.  This can be helpful when determining parameters for survival, crowding
#' and the habitat surface.
#' @param p.seq  Vector sequence of probabilities
#' @param year Ordered Vector sequence of years
#' @param ylim Default y axis scale
#' @param ... Other arguments to \code{plot}
#' @examples plot.survival()
#' plot.survival(xlim=c(0,20),ylim=c(0,0.5))
#'
plot.survival <- function(p.seq=seq(from=1,to=0.5,by=-0.05),years=1:40,ylim=c(0,1),cex=1,...)
{
  plot(x=years,y=p.seq[1]^years,main="Probability of Survival",xlab="Time Steps",type='l',
       ylab="Prob. of survival",ylim=ylim,...)
  text(x=5,y=p.seq[1]^4.5,round(p.seq[1],2))

  for(i in 2:length(p.seq))
  {
    lines(x=years,y=p.seq[i]^years)
    text(x=5,y=p.seq[i]^4.5,round(p.seq[i],2),cex=cex)
  }
}







