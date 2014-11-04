
###########################################################################################
# Some graphics functions

#' Plot the population network in space
#'
#' Plot showing the colonies, and the edges between them, with line thickness indicating weight.
#'   
#'@inheritParams seedPathogen
#'@param lwd Relative line width. This value is scaled to a reasonable value.
#'
#'@name plotColonyNet
#'@family viz 
#'@export

plotColonyNet <- function(pop, lwd=2){
	assert_that(is.numeric(pop$locations), all(c('locations', 'models') %in% names(pop)))
				
	E <- edgeList(pop)
	edgeLocations <- cbind(pop$locations[E[,1],], pop$locations[E[,2],], lwd*E[,3]/max(E[,3]))			

	
	par(mar=c(5,5,1,1)+0.3)
	plot(pop$locations, pch=16, col=brewer.pal(3,'Set1')[1], cex=1,
		ylab='lat', xlab='lon')
	
	apply(edgeLocations,1, function(x) lines(x[c(1,3)],x[c(2,4)], lwd=x[5]))
  points(pop$locations, pch=16, col=brewer.pal(3,'Set1')[1], cex=3)
}


#' A nice heat map for weighted adjacency matrices
#'
#' Plot showing the weight of edges by colouring the adjacency matrix.
#'   
#'@param weightMatrix A matrix 
#'
#'@name plotColonyNet
#'@family viz 

colonyHeat <- function(weightMatrix){
	heatmap(weightMatrix,  scale="none",lwd=3, labRow='', labCol='',  
		col=colorRampPalette(brewer.pal(9,"YlOrRd"))(256))
}

#' A nice heat map for weighted adjacency matrices
#'
#' Plot showing the weight of edges by colouring the adjacency matrix.
#'   
#'@inheritParams seedPathogen
#'
#'@name popHeat
#'@family viz 
#'@export

popHeat <- function(pop){
	colonyHeat(pop$adjacency)
}



###########################################################################################################

# Time series plots


#' Plot the number of susceptible individuals in 

pSus <- function(pop, o = FALSE){
  
  d <- cbind(t(pop$I[1,,]), cumsum(pop$waiting)) %>% data.frame
  colnames(d)[NCOL(d)] <- 'time'

  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longd$value)
  }

  greySelection <- grey(seq(0.2, 0.6, length.out = pop$parameters['nColonies']))

  longd <- melt(d, id = 'time', variable_name = 'Colony')

  ggplot(data = longd,
       aes(x = time, y = value, colour = Colony)) +
    geom_line() +
    ylab('Individuals') + 
    theme_minimal() +
    scale_color_manual(values = greySelection) +
    scale_y_continuous(limits = c(ymin, max(longd$value))) +
    theme(legend.position="none")



}








