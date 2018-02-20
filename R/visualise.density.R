#############################################
# Spatial Genetic Model
#############################################
# Used to visualise the density patterns for a random population over
# a fixed window size.  Useful when setting crowding parameters.
#############
library(spatstat)
library(raster)
###############################################################################
# Set window as 3x3 centred on origin.
# Place the individual we are interested in at (0,0)
# and generate increasing numbers of random individuals in the 3x3 window.
# Store density at (0,0) for specified bandwidth.
###############################################################################

#' Create a table of densities about the origin
#' @description For a given window, randomly place \code{1:num.pts} in a window
#' and measure the density at (0,0) using the bandwidth \code{sigma}.  This can be
#' used to understand the density for randomly distributed populations when
#' setting the crowding parameters for the model.
#' @param win Spatial window of type \code{owin}
#' @param num.pts Maximum number of points to randomly place in the window
#' @param trials Number of times to resample points
#' @param sigma The bandwidth used for \code{density.ppp}
#' @return Matrix with number of rows as points increase from \code{1:num.pts}, each column
#' is the measured density at (0,0) for the particular number of points for each trial.
#' @examples
#' # Create a window centred on (0,0)
#' win <- owin(xrange=c(-2,2),yrange=c(-2,2),unitname="metres")
#' den <- origin.density(win,num.pts=20,trials=50)
#' boxplot(t(den),xlab="Number Points",ylab="Density")

origin.density <- function(win,num.pts=20, trials=10, sigma=1)
{
	xlen <- win$x[2]-win$x[1]
	ylen <- win$y[2]-win$y[1]

	res <- matrix(nrow=num.pts,ncol=trials)
	for(n in 1:num.pts)
	{
		for (trial in 1:trials)
		{
			# Make random points, then add (0,0) as first in list.
			x.pts <- (runif(n)*xlen)-win$x[2]
			y.pts <- (runif(n)*ylen)-win$y[2]
			x.pts <- c(0,x.pts)
			y.pts <- c(0,y.pts)

			p <- ppp(x=x.pts,y=y.pts,window=win)
			d.map <- density.ppp(p,sigma=sigma,at="points",leaveoneout=TRUE)
			res[n,trial] <- d.map[1] # The origin
		}
	}
    res
}
#' Density plot for a given number of points and bandwidth
#' @description This function plots a randomly placed set of \code{num.pts} in
#' the window \code{win}, and colours the surface with the density using the bandwidth
#' \code{sigma}.  The density at (0,0) is displayed.  This is useful when understanding
#' and determining the crowding density and the relationship to bandwidth.
#' @param win Spatial window
#' @param num.pts Number of points randomly placed in the window
#' @param sigma Bandwidth
#' @examples
#' # Create a window centred on (0,0)
#' win <- owin(xrange=c(-2,2),yrange=c(-2,2),unitname="metres")
#' plot.density.example(win,num.pts=10,sigma=1)
#'
plot.density.example <- function(win,num.pts=20,sigma=1)
{
	xlen <- win$x[2]-win$x[1]
	ylen <- win$y[2]-win$y[1]
	x.pts <- (runif(num.pts)*xlen)-win$x[2]
	y.pts <- (runif(num.pts)*ylen)-win$y[2]
	x.pts <- c(0,x.pts)
	y.pts <- c(0,y.pts)

	p <- ppp(x=x.pts,y=y.pts,window=win)
	d.map <- density.ppp(p,sigma=sigma,at="points",leaveoneout=TRUE)
    d.all <- density.ppp(p,sigma=sigma,leaveoneout=TRUE)
	main=paste("Density (0,0): ",round(d.map[1],3),"L:",round(p$n/area.owin(win),3))
	plot(d.all,main="")
	mtext(main,side=3,line=0)
	points(p,pch=16,col='green')
	v.seq = seq(from=win$x[1],to=win$x[2],by=1)
	abline(v=v.seq,lty=2,col='white')
	h.seq = seq(from=win$y[1],to=win$y[2],by=1)
	abline(h=h.seq,lty=2,col='white')
	p.origin <- p[1,]
	points(p.origin,col='white',pch=16,cex=2)
}
# Make the 6x6 window centred on (0,0)
#
win <- owin(xrange=c(-2,2),yrange=c(-2,2),unitname="metres")

# Create table each row is measure of density for particular number of points.
#
#d.1 <- pt.density(win,num.pts=100,trials=10,bw=1)
#d.2 <- pt.density(win,num.pts=100,trials=10,bw=2)
#d.3 <- pt.density(win,num.pts=100,trials=10,bw=3)

