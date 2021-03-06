#############################################
# Spatial Individual-based Model
#############################################
# habitat.r
#############
library(spatstat)
library(raster)
################################################################
# Example operations to define the habitat pixel image
##############################################################

#' Define a circular habitat
#' @description Define a circular habitat using an exponential decay function.
#' The origin(0,0) has value 1, with habitat locations away from the origin
#' decaying in value towards zero.
#'@param x x coordinate
#'@param y y coordinate
#'@param b Decay rate for any habitat location defined by the distance from origin (0,0)
#'@return Habitat value for each location.  This function is used in conjunction with
#'\code{\link[spatstat]{as.im}} to produce an image that defines the habitat.
#'@examples
#' # Example circle habitat for a rectangular window.
#' # hab.win is the habitat window
#' hab.win <- owin(c(-50,50),c(-50,50),unitname="metres")
#' Z <- as.im(circle.habitat,hab.win,b=0.0005)
#' plot(Z)
#'
circle.habitat <- function(x,y,b)
{
	d <- x^2 + y^2
	exp(-b*d)
}
#' Define a set of circular point habitats
#' @description Define a habitat locations within a window using a density
#' function.  The maximum value is 1, and the minimum value in the habitat is
#' set by the argument min.d.
#'@param x x coordinates defining the point locations
#'@param y y coordinates of point locations
#'@param hab.win The habitat window as an owin object
#'@param min.d Minimum value of final habitat map.
#'@param ... Additional arguments to density.ppp
#'@return Habitat value for each location as an im object (\code{\link[spatstat]{as.im}}).
#' This can be used as the habitat window in spatibm.
#'@examples
#' # Example habitat with 3 good habitat locations.
#' # hab.win is the habitat window.  Minimum value is 0.7
#' hab.win <- owin(c(-50,50),c(-50,50),unitname="metres")
#' Use 3 point locations within the -50x50 window
#' x = c(-20,0,20)
#' y = c(-20,0,20)
#' # Note that Z is type im, and sigma,edge are arguments to the density.ppp function.
#' Z <- circle.pts.density(x,y,hab.win,min.d=0.7,sigma=2,edge=FALSE)
#' plot(Z)
#'
circle.pts.density <- function(x,y,hab.win,min.d=0.7,...)
{
  p <- ppp(x=x,y=y,window=hab.win)
  p.d <- density(p,...)
  p.d <- ((p.d-min(p.d))/(max(p.d)-min(p.d)))  # normalise between 0 - 1
  pdm <- as.matrix(p.d)
  pdm[which(pdm < min.d)] <- min.d # set minimum value
  as.im(pdm,W=hab.win)
}
#' Define a circular habitat with constant value
#'@description Create a circular habitat with a value of \code{inside} within a defined distance
#' from the origin(0,0), and \code{outside} elsewhere.
#'@param x x coordinate
#'@param y y coordinate
#'@param dist.inside Distance from the origin (0,0) with habitat value inside
#'@param inside Value of habitat within dist.inside from the origin
#'@param outside Value of habitat further than dist.inside from the origin
#'@return Habitat value for each location.  This function is used in conjunction with
#'\code{\link[spatstat]{as.im}} to produce an image that defines the habitat.
#'@examples
#' hab.win <- owin(c(-50,50),c(-50,50),unitname="metres")
#' Z <- as.im(circle.only,hab.win,dist.inside=10,inside=0.8,outside=0.2)
#' plot(Z)

circle.only <- function(x,y,dist.inside=50,inside=1.0,outside=0.0)
{
	d <- sqrt(x^2 + y^2)
	ifelse(d < dist.inside,inside,outside)
}

######################################################################################
# points.habitat
# IN: pts - a matrix of rows of points representing the centre of good habitat
# If (x,y) is within the threshold.dist of any of the pts then the location is given
# a value of inside, otherwise outside.  Alternatively, could do it as a bandwidth like
# circle habitat...
##################################################################################
loc.dist <- function(x,pts,dist.inside,inside,outside)
{
	pts <- rbind(pts,x)
	dd <- as.matrix(dist(pts))
	diag(dd) <- Inf
	ifelse((min(dd[nrow(dd),]) < dist.inside),inside,outside)
}


#' Define a set of locations with a circular habitat
#' @description Given a set of points (x,y), create a circular habitat centred on
#' each point, with value \code{inside} within \code{dist.inside} and value \code{outside} otherwise.
#' @param x x coordinate
#' @param y y coordinate
#' @param pts Two column matrix or dataframe of x and y coordinates for a set of points
#' @param dist.inside Distance from each point to set value \code{inside}
#' @param inside Value to set for locations within \code{dist.inside} of a point
#' @param outside Value to set for locations further than \code{dist.inside} from a point.
#' @return Habitat value for each location.  This function is used in conjunction with
#'\code{\link[spatstat]{as.im}} to produce an image that defines the habitat.
#'@examples
#' pts <- matrix(data=c(-60,3,-80,-77,30,30,25,10,-40,-30),ncol=2,byrow=T)
#' hab.win <- owin(c(-100,100),c(-100,100),unitname="metres")
#' Z <- as.im(points.habitat,hab.win,pts,dist.inside=20,inside=1.0,outside=0.9)
#' plot(Z)

points.habitat <- function(x,y,pts,dist.inside,inside=1.0,outside=0.8)
{
	# x,y are the vectors of all points
	allxy <- as.matrix(cbind(x,y))
	apply(allxy,1,loc.dist,pts,dist.inside,inside,outside)
}
