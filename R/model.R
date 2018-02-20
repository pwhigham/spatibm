#############################################
# Spatial Individual-based (Genetic) Model
#############################################
# model.R
#############
# Functions to do a single step or multiple steps of a population
#
######################################################################
# update.ids
# IN:curr.pop - current population
#   next.id - next id
# OP: If there are any new borns, update their ids.
# OUT: The new next.id
#######################################################################
update.ids <- function(curr.pop,next.id)
{
	if (curr.pop$n > 0)
	{
		born <- which(curr.pop$marks[,"age"]==1)
		if (length(born) > 0) # we have some kids
		{
			curr.pop$marks[born,"id"] <- as.numeric(next.id:(next.id+length(born)-1))
			next.id <- next.id + length(born)
		}
	}
	list(curr.pop,next.id)
}
#' Perform a single (generation) step for the population
#' @description Perform a single step for a population based on their defined breeding, survival,
#' and crowding parameters and the habitat surface.  This function is only called when the
#' user requires fine-grained control over changes to the population each generation.
#' @param curr.pop The current population as a .ppp
#' @param move.table Movement table as defined by \code{movement.table}.
#' @param survive.table Survival table as defined by \code{survival.table}.
#' @param breed.table Breeding table as defined by \code{breeding.table}
#' @param habitat.surface The habitat (probability) surface. See \code{circle.habitat} for an example, although
#' this surface could be produced in many ways.  The only requirement is that it must cover the entire window
#' associated with the population.
#' @param crowd.table The crowding table as defined by \code{crowding.table}.  May be NA for no crowding.
#' @param crowding.sigma Bandwidth used for density calculation for crowding.  Value of 0.0 indicates no crowding.
#' @param max.dist Maximum distance between breeding parents
#' @return The next population as a .ppp
#'
single.step <- function(curr.pop,
						move.table,
						survive.table,
						breed.table,
						habitat.surface,
						crowd.table=NA, # Default to no crowding
						crowding.sigma=0.0, # Default to no crowding
						max.dist=5.0)
{
	# Who can breed?
	bb <- baseline.breeding(curr.pop,breed.table)
	if (bb$n > 0)  # we have some breeding individuals...
	{
		# match breeding pairs based on closest pairings and maximum separation
		b.pairs <- breeding.pairs(bb,max.dist)
		if (length(b.pairs) > 0)
		{
	# Offspring creation
			children <- breed(bb,b.pairs, breed.table)
			curr.pop <- superimpose(curr.pop,children)  # add children to population
		}
	}
	# Move
	next.pop <- move(curr.pop,move.table)
	# Survival
	next.pop <- baseline.survival(next.pop,survive.table,habitat.surface)
	# and crowding...
	if ((length(crowd.table) > 1) & (crowding.sigma > 0)) next.pop <- crowding.survival(next.pop,crowd.table,sigma=crowding.sigma)
	# Age adjustment
	if (next.pop$n > 0)	next.pop$marks[,"age"] <- next.pop$marks[,"age"]+1 # Update age

    next.pop
}
#' Perform multiple steps (generations) for the population
#' @description Perform \code{steps} generations for a population based on their defined breeding, survival,
#' and crowding parameters and the habitat surface.
#' @param steps Number of steps (generations) to excute the model
#' @param curr.pop The current population as a .ppp
#' @param move.table Movement table as defined by \code{movement.table}.
#' @param survive.table Survival table as defined by \code{survival.table}.
#' @param breed.table Breeding table as defined by \code{breeding.table}
#' @param habitat.surface The habitat (probability) surface. See \code{circle.habitat} for an example, although
#' this surface could be produced in many ways.  The only requirement is that it must cover the entire window
#' associated with the population.  This surface is constant over all generations.  See \code{multiple.step.habitat}
#' for control over a changing habitat for multiple steps.
#' @param crowd.table The crowding table as defined by \code{crowding.table}
#' @param crowding.sigma Bandwidth used for density calculation for crowding
#' @param max.dist Maximum distance between breeding parents
#' @param track.ids Do you want to update the individual ids each generation?  This is only required if you
#' need to produce a family tree.
#' @param trace.output Show a text-based trace of the model running in the console.
#' @return A list containing the population for each generation.  The initial population \code{curr.pop}
#' is the first element in the list.
#'
#'

#############################################################################
# multiple.step
#
# Run model for multiple steps and return the list of populations.
##############################################################################
multiple.step <- function(steps=10,
						curr.pop,
						move.table,
						survive.table,
						breed.table,
						habitat.surface=NA,
						crowd.table=NA,
						crowding.sigma=0.0,  # Default to no crowding
						max.dist=5.0,
						track.ids=FALSE,
						trace.output=TRUE)
{
	pop.series <- vector("list", steps)
	pop.series[[1]] <- curr.pop

	# set next.id from the start population -
	#
	next.id <- max(curr.pop$marks$id) + 1
	for (i in 1:steps)
	{
		if (trace.output)
		{
			cat(paste("Gen:",i," N=",curr.pop$n,"...",sep=""))
			flush.console();
		}
		if (trace.output) if ((i %% 10)==0) cat("\n")

		curr.pop <- single.step(curr.pop,
				 move.table,
				 survive.table,
				 breed.table,
				 habitat.surface,
				 crowd.table,
				 crowding.sigma,
				 max.dist)

		if (track.ids)  # only required if want to build parentage graph
		{
			update.pop.ids <- update.ids(curr.pop,next.id)
			curr.pop <- update.pop.ids[[1]]
			next.id <- update.pop.ids[[2]]
		}

		pop.series[[(i+1)]] <- curr.pop
	}
	if (trace.output) cat("\n")
	pop.series
}

#' Perform multiple steps (generations) for the population, allowing for a stochastic change in habitat.
#' @description Perform \code{steps} generations for a population based on their defined breeding, survival,
#' and crowding parameters and the habitat surface.
#' The \code{habitat.list} defines a list of habitats with a
#' N(m,sd) for each habitat used to determine the duration this habitat is observed.
#' A list of these habitat,N(m,sd) structures defines the temporal pattern of change,
#' and the order of change.  i.e. a list of ( (H1,N1,(H2,N2), ...) defines an ordering of
#' habitat surface.  The N(m,sd) is used to determine the time period before moving from Hi to H(i+1).
#' In addition, there is an implied cycling behavior since when the final (H,N(m,sd)) pattern is
#' used, the next habitat is (H1,N1).
#' @param steps Number of steps (generations) to excute the model
#' @param curr.pop The current population as a .ppp
#' @param move.table Movement table as defined by \code{movement.table}.
#' @param survive.table Survival table as defined by \code{survival.table}.
#' @param breed.table Breeding table as defined by \code{breeding.table}
#' @param habitat.list The list of habitats (probability) surfaces with associated N(m,sd).
#' @param crowd.table The crowding table as defined by \code{crowding.table}
#' @param crowding.sigma Bandwidth used for density calculation for crowding
#' @param max.dist Maximum distance between breeding parents
#' @param track.ids Do you want to update the individual ids each generation?  This is only required if you
#' need to produce a family tree.
#' @param trace.output Show a text-based trace of the model running in the console.
#' @return A list containing the population for each generation.  The initial population \code{curr.pop}
#' is the first element in the list.
#'
multiple.step.habitat <- function(steps=10,
						curr.pop,
						move.table,
						survive.table,
						breed.table,
						habitat.list,
						crowd.table,
						crowding.sigma=1.0,
						max.dist=5.0,
						track.ids=TRUE,
						trace.output=TRUE)
{
	pop.series <- vector("list", steps)
	pop.series[[1]] <- curr.pop

	#
	# Setting up the initial habitat and time before swapping...
	#
	if (class(habitat.list)!="list") stop("Habitat list must be a list")
	if (class(habitat.list[[1]])!="list") stop("Habitat list must be a list of lists")

	curr.habitat.list <- habitat.list[[1]]
	curr.timestep <- 0 # Number of timesteps using current habitat
	curr.habitat <- 1  # Index into the habitat list
	habitat.surface <- curr.habitat.list[[1]]  # The habitat surface
	habitat.maxtime <- ceiling(rnorm(1,curr.habitat.list[[2]],curr.habitat.list[[3]])) # Time

	# set next.id from the start population -
	#
	next.id <- max(curr.pop$marks$id) + 1
	for (i in 1:steps)
	{
		if (trace.output)
		{
			cat(paste("Gen:",i," N=",curr.pop$n,"...",sep=""))
			flush.console();
		}
		if (trace.output) if ((i %% 10)==0) cat("\n")

		curr.pop <- single.step(curr.pop,
				 move.table,
				 survive.table,
				 breed.table,
				 habitat.surface,
				 crowd.table,
				 crowding.sigma,
				 max.dist)

		if (track.ids)  # only required if want to build parentage graph
		{
			update.pop.ids <- update.ids(curr.pop,next.id)
			curr.pop <- update.pop.ids[[1]]
			next.id <- update.pop.ids[[2]]
		}

		pop.series[[(i+1)]] <- curr.pop

		# and deal with habitat surface
		#
		curr.timestep <- curr.timestep + 1
		if (curr.timestep >= habitat.maxtime)  # Time to change...
		{
			curr.habitat <- curr.habitat + 1
			if (curr.habitat > length(habitat.list)) curr.habitat <- 1
			curr.habitat.list <- habitat.list[[curr.habitat]]
			curr.timestep <- 0 # Number of timesteps using current habitat
			habitat.surface <- curr.habitat.list[[1]]  # The habitat surface
			habitat.maxtime <- ceiling(rnorm(1,curr.habitat.list[[2]],curr.habitat.list[[3]])) # Time
		}
	}
	if (trace.output) cat("\n")
	pop.series
}

