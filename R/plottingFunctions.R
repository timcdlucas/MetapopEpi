
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


#' Plot the number of susceptible individuals in each colony
#'
#'@param pop A population object
#'@param o Logical whether to plot from the origin or not
#'@name pSus
#'@export

pSus <- function(pop, start = 1, end = NULL, o = FALSE){
  
  # Just declare these to avoid CRAN check notes
  value <- Colony <- NULL
  
  
  if(is.null(end)){
    end <- pop$parameters['events']/pop$parameters['sample'] + 1
  }

  d <- cbind(t(pop$sample[1, , start:end]), cumsum(pop$sampleWaiting[start:end])) %>% data.frame
  colnames(d)[NCOL(d)] <- 'time'


  greySelection <- grey(seq(0.2, 0.6, length.out = pop$parameters['nColonies']))

  longd <- melt(d, id = 'time', variable.name = 'Colony')
  longd <- longd[longd$value != 0, ]
  
  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longd$value)
  }

  ggplot(data = longd,
       aes(x = time, y = value, colour = Colony)) +
    geom_line() +
    ylab('Individuals') + 
    xlab('Time') + 
    theme_minimal() +
    scale_color_manual(values = greySelection) +
    scale_y_continuous(limits = c(ymin, max(longd$value))) +
    theme(legend.position="none")

}

#' Plot the total number of infected individuals in each colony
#'
#'@inheritParams pClass
#'@name pInf
#'@export


pInf <- function(pop, start = 1, end = NULL, o = FALSE){

  # Just declare these to avoid CRAN check notes
  value <- Colony <- . <- NULL

  if(is.null(end)){
    end <- pop$parameters['events']/pop$parameters['sample'] + 1
  }
  
  I <- pop$sample[2:NROW(pop$sample),,] %>% apply(., c(2,3), sum) %>% t

  d <- cbind(I, cumsum(pop$sampleWaiting)) %>% data.frame
  colnames(d)[NCOL(d)] <- 'time'


  greySelection <- grey(seq(0.2, 0.6, length.out = pop$parameters['nColonies']))

  longd <- melt(d, id = 'time', variable.name = 'Colony')
  longd <- longd[longd$value != 0, ]

  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longd$value)
  }

  ggplot(data = longd,
       aes(x = time, y = value, colour = Colony)) +
    geom_line() +
    ylab('Individuals') + 
    xlab('Time') + 
    theme_minimal() +
    scale_color_manual(values = greySelection) +
    scale_y_continuous(limits = c(ymin, max(longd$value))) +
    theme(legend.position="none")

}







#' Plot the total number of individuals in each colony
#'
#'@inheritParams pClass
#'@name pPop
#'@export


pPop <- function(pop, start = 1, end = NULL, o = FALSE){
  
  # Just declare these to avoid CRAN check notes
  value <- Colony <- . <- NULL

  if(is.null(end)){
    end <- pop$parameters['events']/pop$parameters['sample'] + 1
  }

  I <- pop$sample[, , start:end] %>% apply(., c(2,3), sum) %>% t

  d <- cbind(I, cumsum(pop$sampleWaiting[start:end])) %>% data.frame
  colnames(d)[NCOL(d)] <- 'time'


  greySelection <- grey(seq(0.2, 0.6, length.out = pop$parameters['nColonies']))

  longd <- melt(d, id = 'time', variable.name = 'Colony')
  longd <- longd[longd$value != 0, ]

  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longd$value)
  }

  ggplot(data = longd,
       aes(x = time, y = value, colour = Colony)) +
    geom_line() +
    ylab('Individuals') + 
    xlab('Time') + 
    theme_minimal() +
    scale_color_manual(values = greySelection) +
    scale_y_continuous(limits = c(ymin, max(longd$value))) +
    theme(legend.position="none")

}



#' Plot something...
#'
#'@inheritParams pClass
#'@name pAll
#'@export


pAll <- function(pop, start = 1, end = NULL, o = FALSE){

  # Just declare these to avoid CRAN check notes
  value <- colonyState <- . <- NULL

  if(is.null(end)){
    end <- pop$parameters['events']/pop$parameters['sample'] + 1
  }

  # Sum infections from different classes
  z <- lapply(1:pop$parameters['nPathogens'], function(x) apply(pop$sample[pop$whichClasses[, x], , start:end], c(2, 3), sum))

  z2 <- do.call(rbind, z)

  
  d <- cbind(t(z2), cumsum(pop$sampleWaiting[start:end])) %>% data.frame
  d <- cbind(d, 'disease')
  colnames(d)[c(NCOL(d) - 1, NCOL(d))] <- c('time', 'state')

  s <- cbind(t(pop$sample[1, , start:end]), cumsum(pop$sampleWaiting[start:end])) %>% data.frame
  s <- cbind(s, 'susceptible')
  colnames(s)[c(NCOL(s) - 1, NCOL(s))] <- c('time', 'state')

  longd <- melt(d, id = c('time', 'state') , variable.name = 'colony')
  longs <- melt(s, id = c('time', 'state') , variable.name = 'colony')
  
  longAll <- rbind(longd, longs)

  longAll$colonyState <- paste0(longAll$colony, longAll$state) %>% factor
  longAll <- longAll[longAll$value != 0, ]

  twoGroups <- rep(brewer.pal(3, 'Set1')[2], length(table(longAll$colonyState)))
  twoGroups[grep('disease', names(table(longAll$colonyState)))] <- brewer.pal(3, 'Set1')[1]
  names(twoGroups) <- names(table(longAll$colonyState))

  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longAll$value)
  }

  ggplot(data = longAll,
       aes(x = time, y = value, colour = colonyState)) +
    geom_line() +
    ylab('Individuals') + 
    xlab('Time') + 
    theme_minimal() +
    scale_color_manual(values = twoGroups) +
    scale_y_continuous(limits = c(ymin, max(longAll$value))) +
    theme(legend.position="none")

}


#' Plot the total number of individuals in each disease class
#'
#'@param pop A MetapopEpi class object
#'@param start What event to start plotting from.
#'@param end What event to end plotting from. If NULL, plot until end.
#'@param S Logical, controls whether or not to include the susceptible population.
#'@param nPath Logical, Colour the plot by number of pathogens or by class number
#'@param o Logical. Truncate y-axis.
#'@name pClass
#'@export

pClass <- function(pop, start = 1, end = NULL, S = TRUE, nPath = TRUE, o = FALSE){
  
  # Just declare these to avoid CRAN check notes
  value <- colony <- disease <- NULL

  if(is.null(end)){
    end <- pop$parameters['events']/pop$parameters['sample'] + 1
  }

  if(S){
    removeS <- 1:pop$nClasses
    nClass <- pop$nClasses
  } else {
    removeS <- -1
    nClass <- pop$nClasses - 1
  }

  # Sum infections from different classes
  z <- t(apply(pop$sample[removeS, , start:end], c(1, 3), sum))
   
  
  d <- cbind(z, cumsum(pop$sampleWaiting[start:end])) %>% data.frame

  nPaths <- sapply(pop$diseaseList, length)[removeS]

  colnames(d) <- c(1:(nClass), 'time')
  


  longd <- melt(d, id = 'time' , variable.name = 'disease')

  longd$nPath <- nPaths[longd$disease]

  longd <- longd[longd$value != 0, ]
  
  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longd$value)
  }

  cols <- brewer.pal(8, 'Dark2')

  ggplot(data = longd,
       aes(x = time, y = value, group = disease, colour = factor(nPath))) +
    geom_line() +
    ylab('Individuals') + 
    xlab('Time') + 
    theme_minimal() +
    scale_y_continuous(limits = c(ymin, max(longd$value))) +
    labs(colour = "Coinfection lvl") +
    scale_color_brewer(palette="Dark2")

}





#' Plot the total number of infected individuals with each disease
#'
#'@inheritParams pClass
#'@name pDis
#'@export


pDis <- function(pop, start = 1, end = NULL, o = FALSE){
  
  # Just declare these to avoid CRAN check notes
  value <- colony <- Pathogen <- NULL


  # If null set end to final event (i.e. plot whole time)
  if(is.null(end)){
    end <- dim(pop$sample)[3]
  }


  
  I <- lapply(1:pop$parameters['nPathogens'], function(p) colSums(colSums(pop$sample[pop$whichClasses[, p], , start:end])))

  I <- do.call(cbind, I)

  d <- cbind(I, cumsum(pop$sampleWaiting[start:end])) %>% data.frame
  colnames(d)[NCOL(d)] <- 'time'


  greySelection <- grey(seq(0.2, 0.6, length.out = pop$parameters['nPathogens']))

  longd <- melt(d, id = 'time', variable.name = 'Pathogen')
  longd <- longd[longd$value != 0, ]


  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longd$value)
  }

  ggplot(data = longd,
       aes(x = time, y = value, colour = Pathogen)) +
    geom_line() +
    ylab('Individuals') + 
    xlab('Time') + 
    theme_minimal() +
    scale_color_manual(values = greySelection) +
    scale_y_continuous(limits = c(ymin, max(longd$value))) +
    theme(legend.position="none")

}



#' Plot the total number of infected individuals and total Susceptible.
#'
#'@inheritParams pClass
#'@name pSI
#'@export

pSI <- function(pop, start = 1, end = NULL){
  
  # Just declare these to avoid CRAN check notes
  value <- SI_class <- NULL

  if(is.null(end)){
    end <- pop$parameters['events']/pop$parameters['sample'] + 1
  }

  I <- cbind(apply(pop$sample[-1, , start:end], 3, sum), apply(pop$sample[1, , start:end], 2, sum), cumsum(pop$sampleWaiting[start:end])) %>% data.frame
  colnames(I) <- c('I', 'S', 'time')

  longd <- melt(I, id = 'time', variable.name = 'SI_class')
  longd <- longd[longd$value != 0, ]

  ggplot(data = longd,
       aes(x = time, y = value, colour = SI_class)) +
    geom_line() +
    ylab('Individuals') + 
    xlab('Time') + 
    theme_minimal() +
    scale_y_continuous(limits = c(0, max(longd$value))) +
    theme(legend.position="none")

}


