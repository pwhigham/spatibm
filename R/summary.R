#############################################
# Spatial Individual-based (Genetic) Model
#############################################
# summary.R
#############
####################################################################
# Example summary functions
####################################################################

# This function is internal
#
count.pop <- function(p,x.col,x.val)
{
	if (is.null(p)) return(0)
	if (x.col==0) return(p$n)

	length(which(p$marks[,x.col]==x.val))
}
# Plot total pop, M/F
#
#' Plot total population
#' @description Create a single plot showing total number of individuals,
#' and number of males/females in the population for each generation.
#' @param p.list  List of populations created by running a model.
#
plot.pop.count <- function(p.list)
{
	total.n <- unlist(lapply(p.list,count.pop,x.col=0,x.val=0))
	m.n <- unlist(lapply(p.list,count.pop,x.col=2,x.val=1))
	f.n <- unlist(lapply(p.list,count.pop,x.col=2,x.val=2))

	plot(x=1:length(p.list), y=total.n,type='l', lwd=1.5,
			main="",xlab="Time Step",ylab="Pop Count",
			ylim=c(min((m.n),(f.n)),max(total.n)))

	points(x=1:length(p.list),y=total.n,pch=16,col=1,cex=0.7)

	lines(x=1:length(p.list),m.n,lty=1,col=2,lwd=0.75)
	lines(x=1:length(p.list),f.n,lty=2,col=3,lwd=0.75)

	points(x=1:length(p.list),m.n,pch=16,col=2,cex=0.7)
	points(x=1:length(p.list),f.n,pch=16,col=3,cex=0.7)

	legend("topleft",c("All","M","F"),col=1:3,pch=c(16,16,16),cex=0.7)
}
#
# Calculate Wright F statistic
# for each loci.
#' Calculate Wright F statistic for a population
#' @param p Population
#' @return Numeric vector with 4 values per loci.  Values are p frequency, H observed,
#' H expected and the F statistic.
#' @examples
#' p <- create_ibm_population()
#' # The default population has two loci, so there are 2 x 4 values.
#' print(F.stat(p))
#'
F.stat <- function(p)
{
  first <- which(colnames(p$marks)=="V1") # First of the 2nd columns
	fstats <- NULL

	while(first < ncol(p$marks))
	{
		alleles <- p$marks[,first:(first+1)]
		numAA <- length(which(apply(alleles,1,sum)==2)) #AA
		numaa <- length(which(apply(alleles,1,sum)==0)) #aa
		numAa <- length(which(apply(alleles,1,sum)==1)) # Aa

		p.freq <- ((2*numAA)+numAa)/(2*p$n)
		q.freq <- 1 - p.freq
		Hobs <- numAa/p$n
		Hexp <- 2*p.freq*q.freq
		Fs <- (Hexp-Hobs)/Hexp
		fstats <- c(fstats,p.freq,Hobs,Hexp,Fs)
		first <- first + 2
	}
	fstats
}
#' Population summary statistics
#' @description Given a single population as .ppp, construct a vector
#' that summarises the properties of the population.
#' @param pop The population as a .ppp
#' @param gen The generation of this population.
#' @param sigma The bandwidth to use for calculating point density of individuals
#' @return Numeric vector representing a summary of population properties.
#' See \code{model.summary} to create a labelled table with this summary data.
pop.summary <- function(pop,gen,sigma=1.0)
{
	# Num-inds, mean age, youngest, oldest, malesexratio, mean-density, min-den, max-den
	#
	if (pop$n == 0) return (c(gen,pop$n,rep(0,8)))

	den <- density.ppp(pop,sigma=sigma,at="points",leaveoneout=TRUE)

	c(gen,pop$n,mean(pop$marks$age),min(pop$marks$age),max(pop$marks$age),
			(length(which(pop$marks$sex==1))/pop$n),mean(den),min(den),max(den),
			F.stat(pop))
}
#' Create summary table for a model run
#' @description Given a list of populations, probably produced from running
#' a model (see \code{multiple.step}), create a summary table, one row per generation.
#' @param p.list A list of point populations (.ppp)
#' @param sigma Bandwidth used for calculating density for each generation
#' @file.csv A file name as a string or NA.  If not NA, then a csv file using
#' this name is created.
#' @return A dataframe, one row per generation, showing summary statistics for the population.
#' @examples
#' # Create a default population
#' p <- create_ibm_population()
#' # Create bigger window for population to expand
#' p$window <- owin(c(-50,50),c(-50,50),unitname="metres")
#' # Use default parameters for behaviour.
#' pm <- multiple.step(steps=5,curr.pop=p,move.table=movement.table(),
#' survival.table(),breeding.table(),habitat.surface=NA,crowding.table())
#' print(model.summary(pm))
model.summary <- function(p.list,sigma=1,file.csv=NA)
{
	res <- NULL

	for (i in 1:length(p.list))
	{
		if (p.list[[i]]$n==0) break;
		res <- rbind(res,pop.summary(p.list[[i]],(i-1),sigma))
	}
	res <- as.data.frame(res)
	# calculate column names for gene frequence - 3 per loci
    first <- which(colnames(p.list[[1]]$marks)=="V1") # First of the 2nd columns
	num.alleles <-((ncol(p.list[[1]]$marks)-first)+1)/2

	# Make string names for allele and Fstat summary data
	fstat.strings <- paste0(c("pfreq","Hobs","Hexp","Fstat"),
	                        sort(rep(1:num.alleles,4))
	                        )
	colnames(res) <- c("time","n","mean.age","min.age","max.age","sexratio",
						"mean.den","min.den","max.den",
						fstat.strings)

	if (!is.na(file.csv)) write.csv(file.csv,res,header=T)
	res
}
#' Plot all generations of a model as points
#' @description Plot all generations of a model run, one plotting window
#' per generation.  Note that the number of plotting windows needs to be
#' set prior to calling (probably using \code{par(mfrow=c(r,c),...}).
#' @param p.list List of populations as .ppp, probably produced by running a model.
#' @param xlim,ylim Required explicitly if used, since \code{plot.ppp} uses the .ppp window
#' as the default window to scale the plot.
#' @param cex Size of points used for plotting
#' @param pch Point type
#' @param ...  Other graphical parameters
#' @examples
#' # Create a default population
#' p <- create_ibm_population()
#' # Create bigger window for population to expand
#' p$window <- owin(c(-50,50),c(-50,50),unitname="metres")
#' # Use default parameters for behaviour.
#' pm <- multiple.step(steps=5,curr.pop=p,move.table=movement.table(),
#' survival.table(),breeding.table(),habitat.surface=NA,crowding.table())
#' op <- par(mfrow=c(2,3),mar=c(1,1,0.5,0.5))  # setup graphics
#' plot.all.gens(pm)
#' par(op) # reset graphics
plot.all.gens <- function(p.list,xlim=NA,ylim=NA,cex=0.7,pch=16,...)
{
	for (i in 1:length(p.list))
	{
		if (length(p.list[[i]])==0) break;
		if (p.list[[i]]$n==0) break;
		pwin <- p.list[[i]]$window

		if (length(xlim) > 1)
		{
			p.list[[i]]$window[[2]] <- xlim
			p.list[[i]]$window[[3]] <- ylim
		}
		plot(p.list[[i]],which.marks=2, # male/female
			clipwin=p.list[[i]]$window,
			cols=c("red","green"),
			cex=cex,pch=pch,
			legend=FALSE,main="",
			 ...)
		mtext(i,side=3,line=-2)
		p.list[[i]]$window <- pwin
	}
}

