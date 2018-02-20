#####################################################
# Use igraph to draw family tree of model results
#####################################################

# Internal function
single.mark.table <- function(sp)
{
	res <- NULL
	for (i in 1:length(sp))
	{
		res <- rbind(res,sp[[i]]$marks)
	}
	res
}
#' Create a tree-style layout for igraph
#' @description Creates a layout design for a graph that assumes we
#' are showing a family tree structure.  Creates a layout for boxes to be
#' drawn using the \code{allelebox} object for each individual.
#' @param sp Graph created using \code{build.population.graph}
#' @param w Width of allelebox
#' @param h Height of allelebox
#' @param w.gap Horizontal gap between boxes in graphical units
#' @param h.gap Vertical gap between boxes in graphical units
#' @return The igraph layout

layout.as.tree <- function(sp,w=1,h=1,w.gap=1,h.gap=1)
{
	res <- NULL
	num.levels <- length(sp)
	start.y <- (num.levels)*(h.gap+h)

	# We are centering each level around zero

	level1 <- sp[[1]]
	if (level1$n==0) return(NULL)
	start.x <- -floor(level1$n/2)*(w+w.gap)

	for (i in 1:level1$n)
	{
		res <- rbind(res,c(start.x,start.y))
		start.x <- start.x + (w+w.gap)
	}
	for (i in 2:length(sp))
	{
		start.y <- start.y - (h.gap+h)
		if (length(sp[[i]]) > 0)
		{
			if (sp[[i]]$n==0) break;

			kids <- length(which(sp[[i]]$marks$age==1))
			if (kids > 0)
			{
				start.x <- -floor(kids/2)*(w+w.gap)
				for (k in 1:kids)
				{
					res <- rbind(res,c(start.x,start.y))
					start.x <- start.x + (w+w.gap)
				}
			}
		}
	}
	res
}
######################################################
# Drawing object for the "allele" type of object
# Used for igraph plotting of family tree.
# Not called by user.
######################################################
allelebox <- function(coords, v=NULL, params)
{
  box.width <- params("vertex", "width")
  bw <- box.width[1]
  box.height <- params("vertex", "height")
  bh <- box.height[1]
  allele.vals <- params("vertex","alleles")
  num.alleles <- length(allele.vals[[1]])

  mapply(coords[,1],coords[,2],allele.vals,
		 FUN=function(x,y,aval){

			a.w <- (bw*2)/num.alleles
			a.h <- bh
			s.w <- x-bw
			for (aspot in 1:num.alleles)
			{
				 polygon(
				 x=c(s.w,s.w+a.w,
					 s.w+a.w,s.w),
				 y=c(y-bh,y-bh,
					 y+bh,y+bh),
				 col=ifelse(aval[aspot]==1,"white","black"))
				 s.w <- s.w + a.w
			 }

		 })
}
#' Build a graph for a list of populations
#' @description Create a graph using the results of a model run.
#' Assumes that p.list is a list of populations as .ppp, probably
#' created using \code{multiple.step}.
#' @param p.list Population list
#' @return A graph as an igraph object
#'
build.population.graph <- function(p.list)
{
	marks <- single.mark.table(p.list)
	u.m <- unique(marks$id)
	g <- make_empty_graph(n=length(u.m))

	kids <- marks[which(marks$age==1),]
	if (nrow(kids) > 0)
	{
		for (k in 1:nrow(kids)) # for each unique kid
		{
		# add edges from parents to kid
			g <- add_edges(g,c(kids[k,"f"],kids[k,"id"]))
			g <- set_edge_attr(g,"color",index=E(g)[length(E(g))],value="red")
			g <- add_edges(g,c(kids[k,"m"],kids[k,"id"]))
			g <- set_edge_attr(g,"color",index=E(g)[length(E(g))],value="black")
		}
	}
	g
}
get.allele.list <- function(sp)
{
	marks <- single.mark.table(sp)
	u.m <- unique(marks$id)
	# assumes that the ids are 1 - N where N is the number of unique
	# individuals.  Need to check next.id allocation...
	result <- vector("list", length(u.m))
	for (i in 1:length(u.m))
	{
		pid <- marks[which(marks$id==u.m[i]),]
		result[[i]] <- as.numeric(pid[1,6:ncol(marks)])
	}
	result
}


# add_shape("alleles", clip=shapes("circle")$clip,
             # plot=allelebox,
			 # parameters=list(vertex.width=1,
							 # vertex.height=1,
							 # vertex.alleles=length(V(g))))

#' Plot the population as a family tree
#' @description Plot the population as a family tree.  Note that prior
#' to using this the alleles shape needs to be added.  See the example below. Edges
#' are red from females, black from males. This function is in development, and only
#' works for small populations and quite a lot of massaging with the parameters and layouts.
#' @param g igraph created using \code{build.population.graph}.
#' @param pop.series The list of populations used to build the population graph.
#' @param tree.layout The tree layout used for plotting
#' @param vertex.height Vertex box height
#' @param vertex.width Vertex box width
#' @param edge.width Edge width of connecting lines
#' @param edge.arrow.size Arrow size of graph
#' @param edge.curved Are the drawn edges curved?
#' @param ... other parameters to plot.igraph
#' @examples
#' # Create a small population of size 10
#' N <- 10
#' p <- create_ibm_population(N)
#' # Create bigger window for population to expand
#' p$window <- owin(c(-50,50),c(-50,50),unitname="metres")
#' # Use default parameters for behaviour.
#' # Note that we require the ids for a tree, so track.ids=TRUE
#' pop.series <- multiple.step(steps=5,curr.pop=p,move.table=movement.table(),
#' survival.table(),breeding.table(),habitat.surface=NA,crowding.table(),
#' track.ids=TRUE)
#' g <- build.population.graph(pop.series)
#'
#' # And now add the shape function - depends on the number of
#' # alleles, so needs to be done after the graph has been built.
#'
#' add_shape("alleles", clip=shapes("circle")$clip,
#'        plot=allelebox,parameters=list(vertex.width=1,vertex.height=1,vertex.alleles=length(V(g))))
#' plot.population.graph(g,pop.series,tree.layout=layout_with_fr(g),asp=1,vertex.height=0.02,vertex.width=0.04)
#
#' # Not run - using spatibm package layout
#' #
#' # plot.population.graph(g,pop.series,tree.layout=layout_as_tree(g),asp=1,vertex.height=0.02,vertex.width=0.04)
#' #

plot.population.graph <- function(g,pop.series,
								tree.layout=NA,
								vertex.height=0.02,
								vertex.width=0.04,
								edge.width=rep(0.02,length(E(g))),
								edge.arrow.size=0.35,
								edge.curved=TRUE,
								...)
{

	gal <- get.allele.list(pop.series)
	if (length(tree.layout)==1) tree.layout <- layout.as.tree(pop.series,
					w=1,h=1,w.gap=2,h.gap=2)


	plot(g,vertex.shape="alleles",vertex.label=NA,
		vertex.height=vertex.height,
		vertex.width=vertex.width,
	   vertex.alleles=gal,
	   edge.width=edge.width,
	   edge.arrow.size=edge.arrow.size,
	   layout=tree.layout,
	   edge.curved=edge.curved, ...)
}
# NOT RUN
#
#g <- build.population.graph(pop.series)

# And now add the shape function - depends on the number of
# alleles, so needs to be done after g has been defined.
#

# add_shape("alleles", clip=shapes("circle")$clip,
             # plot=allelebox,
			 # parameters=list(vertex.width=1,
							 # vertex.height=1,
							 # vertex.alleles=length(V(g))))

#plot.population.graph(g,pop.series,tree.layout=layout_with_fr(g),asp=1,
#				vertex.height=0.02,vertex.width=0.04)
#plot.population.graph(g,pop.series,tree.layout=NA,asp=0,vertex.height=0.02,
#								vertex.width=0.04)

# gal <- get.allele.list(pop.series)
# tree.layout <- layout.as.tree(pop.series,
					# w=1,h=1,w.gap=2,h.gap=2)


# NOTES:
# layout.as.tree - works ok for small networks -- my version so takes the population series.
# layout_with_fr(g) is interesting because it shows how local clustering of genes occurs
# layout_as_tree(g,root=1:N) should have worked ok, but I couldn't stop the vertex drawings from
# overlapping.

# plot(g,vertex.shape="alleles",vertex.label=NA,
	   # vertex.height=0.04,
	   # vertex.width=0.04,
	   # vertex.alleles=gal,
	   # edge.width=rep(0.02,length(E(g))),
	   # edge.arrow.size=0.35,
	   # #layout=layout_as_tree(g,root=1:N),
	   # #layout=layout_in_circle(g),
	   # #layout=layout_with_fr(g),
	   # layout=tree.layout,
	   # #asp=0,
	   # #axes=TRUE,
	   # margin=-0.1,
	   # edge.curved=TRUE)

