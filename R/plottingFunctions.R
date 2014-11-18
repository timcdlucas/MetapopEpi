
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

pSus <- function(pop, o = FALSE){
  
  d <- cbind(t(pop$I[1,,]), cumsum(pop$waiting)) %>% data.frame
  colnames(d)[NCOL(d)] <- 'time'


  greySelection <- grey(seq(0.2, 0.6, length.out = pop$parameters['nColonies']))

  longd <- melt(d, id = 'time', variable.name = 'Colony')

  
  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longd$value)
  }

  ggplot(data = longd,
       aes(x = time, y = value, colour = Colony)) +
    geom_line() +
    ylab('Individuals') + 
    theme_minimal() +
    scale_color_manual(values = greySelection) +
    scale_y_continuous(limits = c(ymin, max(longd$value))) +
    theme(legend.position="none")

}

#' Plot the total number of infected individuals in each colony
#'
#'@inheritParams pSus
#'@name pInf
#'@export


pInf <- function(pop, o = FALSE){
  
  I <- pop$I[2:NROW(pop$I),,] %>% apply(., c(2,3), sum) %>% t

  d <- cbind(I, cumsum(pop$waiting)) %>% data.frame
  colnames(d)[NCOL(d)] <- 'time'


  greySelection <- grey(seq(0.2, 0.6, length.out = pop$parameters['nColonies']))

  longd <- melt(d, id = 'time', variable.name = 'Colony')

  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longd$value)
  }

  ggplot(data = longd,
       aes(x = time, y = value, colour = Colony)) +
    geom_line() +
    ylab('Individuals') + 
    theme_minimal() +
    scale_color_manual(values = greySelection) +
    scale_y_continuous(limits = c(ymin, max(longd$value))) +
    theme(legend.position="none")

}







#' Plot the total number of individuals in each colony
#'
#'@inheritParams pSus
#'@name pPop
#'@export


pPop <- function(pop, o = FALSE){
  
  I <- pop$I %>% apply(., c(2,3), sum) %>% t

  d <- cbind(I, cumsum(pop$waiting)) %>% data.frame
  colnames(d)[NCOL(d)] <- 'time'


  greySelection <- grey(seq(0.2, 0.6, length.out = pop$parameters['nColonies']))

  longd <- melt(d, id = 'time', variable.name = 'Colony')


  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longd$value)
  }

  ggplot(data = longd,
       aes(x = time, y = value, colour = Colony)) +
    geom_line() +
    ylab('Individuals') + 
    theme_minimal() +
    scale_color_manual(values = greySelection) +
    scale_y_continuous(limits = c(ymin, max(longd$value))) +
    theme(legend.position="none")

}



#' Plot something...
#'
#'@inheritParams pSus
#'@name pAll
#'@export


pAll <- function(pop, o = FALSE){

  # Sum infections from different classes
  z <- lapply(1:pop$parameters['nPathogens'], function(x) apply(pop$I[pop$whichClasses[,x],,], c(2,3) , sum))

  z2 <- do.call(rbind, z)

  
  d <- cbind(t(z2), cumsum(pop$waiting)) %>% data.frame
  d <- cbind(d, 'disease')
  colnames(d)[c(NCOL(d)-1, NCOL(d))] <- c('time', 'state')

  s <- cbind(t(pop$I[1,,]), cumsum(pop$waiting)) %>% data.frame
  s <- cbind(s, 'susceptible')
  colnames(s)[c(NCOL(s)-1, NCOL(s))] <- c('time', 'state')

  longd <- melt(d, id = c('time', 'state') , variable.name = 'colony')
  longs <- melt(s, id = c('time', 'state') , variable.name = 'colony')
  
  longAll <- rbind(longd, longs)

  longAll$colonyState <- paste0(longAll$colony, longAll$state) %>% factor


  twoGroups <- rep(brewer.pal(3, 'Set1')[2], length(table(longAll$colonyState)))
  twoGroups[grep('disease',names(table(longAll$colonyState)))] <- brewer.pal(3, 'Set1')[1]
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
    theme_minimal() +
    scale_color_manual(values = twoGroups) +
    scale_y_continuous(limits = c(ymin, max(longAll$value))) +
    theme(legend.position="none")

}


#' Plot the total number of individuals in each disease class
#'
#'@inheritParams pSus
#'@name pClass
#'@export

pClass <- function(pop, o = FALSE){

  # Sum infections from different classes
  z <- lapply(1:pop$parameters['nPathogens'], function(x) colSums(apply(pop$I[pop$whichClasses[,x],,], c(2,3) , sum)))

  z2 <- do.call(cbind, z)

  
  d <- cbind(z2, cumsum(pop$waiting)) %>% data.frame
  #d <- cbind(d, 'disease')
  colnames(d) <- c(paste0('p', 1:pop$parameters['nPathogens']),'time')

  s <- data.frame(time = cumsum(pop$waiting), disease = 's', value = colSums(pop$I[1,,]))

  longd <- melt(d, id = 'time' , variable.name = 'disease')
  
  
  longAll <- rbind(longd, s)

  #longAll$colonyState <- paste0(longAll$colony, longAll$state) %>% factor


  #twoGroups <- rep(brewer.pal(3, 'Set1')[2], length(table(longAll$colonyState)))
  #twoGroups[grep('disease',names(table(longAll$colonyState)))] <- brewer.pal(3, 'Set1')[1]
  #names(twoGroups) <- names(table(longAll$colonyState))

  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longAll$value)
  }

  ggplot(data = longAll,
       aes(x = time, y = value, colour = disease)) +
    geom_line() +
    ylab('Individuals') + 
    theme_minimal() +
    scale_y_continuous(limits = c(ymin, max(longAll$value))) +
    theme(legend.position="none")
  

}





#' Plot the total number of infected individuals with each disease
#'
#'@inheritParams pSus
#'@name pDis
#'@export


pDis <- function(pop, o = FALSE){
  
  I <- lapply(1:pop$parameters['nPathogens'], function(p) colSums(colSums(pop$I[pop$whichClasses[,p], , ])))

  I <- do.call(cbind, I)

  d <- cbind(I, cumsum(pop$waiting)) %>% data.frame
  colnames(d)[NCOL(d)] <- 'time'


  greySelection <- grey(seq(0.2, 0.6, length.out = pop$parameters['nPathogens']))

  longd <- melt(d, id = 'time', variable.name = 'Pathogen')

  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longd$value)
  }

  ggplot(data = longd,
       aes(x = time, y = value, colour = Pathogen)) +
    geom_line() +
    ylab('Individuals') + 
    theme_minimal() +
    scale_color_manual(values = greySelection) +
    scale_y_continuous(limits = c(ymin, max(longd$value))) +
    theme(legend.position="none")

}



#' Plot the total number of infected individuals and total Susceptible.
#'
#'@inheritParams pSus
#'@name pSI
#'@export

pSI <- function(pop, o = FALSE){
  
  I <- cbind(apply(pop$I[-1, , ], 3, sum), apply(pop$I[1, , ], 2, sum), cumsum(pop$waiting)) %>% data.frame
  colnames(I) <- c('I', 'S', 'time')

  if(o){
    ymin <- 0
  } else { 
    ymin <- min(longd$value)
  }

  longd <- melt(I, id = 'time', variable.name = 'SI_class')

  ggplot(data = longd,
       aes(x = time, y = value, colour = SI_class)) +
    geom_line() +
    ylab('Individuals') + 
    theme_minimal() +
    scale_y_continuous(limits = c(ymin, max(longd$value))) +
    theme(legend.position="none")

}


