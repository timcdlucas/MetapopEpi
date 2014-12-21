



#'A package containing functions to run multi pathogen, metapopulation epidemiological simulations.
#'
#'@name MetapopEpi
#'@author Tim CD Lucas
#'@docType package
#'@import assertthat RColorBrewer magrittr igraph ggplot2 reshape2

NULL




#' Add single individuals infected with a number of pathogens
#' 
#' Allows addition of an infected individual of 1 or more pathogens. 
#'
#'@param pop A metapopulation object as created by makePop().
#'@param pathogens Character vector of which pathogens should be seeded
#'@param n How many individuals should be infected.
#'@param diffCols Logical. If True, pathogen infections are forced to be in different colonies. Allows n to be as big as meanColonySize but requires nColonies > nPathogens
#'
#'@return An update object of same structure as from makePop
#'@family Run,sims
#'
#'@export
#'@name seedPathogen
#'@examples
#'p <- makePop()
#'u <- seedPathogen(p, pathogens = 1:p$parameters['nPathogens'])


seedPathogen <- function(pop, pathogens = 1, n = 1, diffCols = TRUE){
	
  assert_that(all(sapply(pathogens, is.count)))
	# can't seed more species than were initiliased.
	assert_that(max(pathogens) <= pop$parameters['nPathogens'])
  assert_that(length(pathogens) <= pop$parameters['nPathogens'])
	
  # choose random colony
	r <- sample(pop$parameters['nColonies'], length(pathogens), replace = !diffCols)  

	for(path in 1:length(pathogens)){

    
		pop$I[pathogens[path] + 1, r[path], 1] <- pop$I[pathogens[path] + 1, r[path], 1] + n
		# keep pop size constant
		pop$I[1, r[path], 1] <- pop$I[1, r[path], 1] - n
	}

  pop$sample[, , 1] <- pop$I[, , 1]
  pop <- transRates(pop, 1)

  assert_that(all(pop$I[, , 1] >= 0))

	return(pop)
}


#' Perform a random event
#' 
#' Randomly selects an event weighted by their rates, runs the event and then calculates a waiting time.
#'@param pop A population object
#'@param t Time step. Note the population should be AT time t, going to t+1.
#'@name randEvent 

randEvent <- function(pop, t, tMod){
  # Calculate random number each time. Expect, precalcing randUnifs might be quicker.


  # Using precalcuated random uniforms in [0,1], pop$randEventU
  # Find which interval this falls in as quick way to select value
  # Should have a think about open and closed intervalsfromClass toColony toClass rate

  event <- pop$transitions[findInterval(pop$randEventU[t] * pop$totalRate, cumsum(c(0, pop$transitions$rate))), ]


  # Copy pop to t + 1
  pop$I[, , tMod + 1] <- pop$I[, , tMod]
  
  # run event
  pop$I[event[['fromClass']], event[['fromColony']], tMod + 1] <- 
    pop$I[event[['fromClass']], event[['fromColony']], tMod ] - 1
  pop$I[event[['toClass']], event[['toColony']], tMod + 1] <- 
    pop$I[event[['toClass']], event[['toColony']], tMod] + 1

  pop <- transRates(pop, tMod + 1)
  
  pop <- waitingTime(pop, tMod)

  return(pop)
}


  

#' Run simulation up to certain event number
#'
#' Each time step, pick a random event and perform that event.

#'@export
#'@name runSim
#'@inheritParams randEvent
#'@param time The event number to run the simulation until.

runSim <- function(pop, time = 'end'){

  assert_that(is.count(time) | time == 'end')

  if (time == 'end'){
    time <- pop$parameters['events']
  }
  

  pb <- txtProgressBar(1, pop$parameters['events'], style = 3)
    
  for (t in 1:time){

    

    tMod <- (t - 1) %% pop$parameters['sample'] + 1


    pop <- randEvent(pop, t, tMod)

    if(tMod == pop$parameters['sample']){
      pop$sampleWaiting[t/pop$parameters['sample'] + 1] <- sum(pop$waiting)
      pop$sample[, , t/pop$parameters['sample'] + 1] <- pop$I[, , tMod + 1]
      pop$I[, , 1] <- pop$I[, , tMod + 1]
    }
  
    setTxtProgressBar(pb, t)
  }

  close(pb)
  return(pop)
}

###############################################################################



#' Takes a number [0,1] (which is presumably random) and returns waiting time.
#'
#' A single draw from a exponential distribution given randU, a number between 0 and 1.
#'   As population now has 
#'
#'@name waitingTime
#'@inheritParams randEvent

waitingTime <- function(pop, t){
  pop$waiting[t + 1] <- log(1 - pop$randU[t]) / (-pop$totalRate)
  return(pop)
}






#' Calculate new transition rates.
#'
#'
#'@inheritParams randEvent
#'@name transRates
#'@family initialRates
#'@export 


transRates <- function(pop, t){

  # Infection from which class to which class.

  rate <- c(birthR(pop, t), deathR(pop, t),  infectionR(pop, t), coinfectionR(pop, t), dispersalR(pop, t))

  if(pop$models$model == 'SIS'){
    rate <- c(rate, recoveryR(pop, t))
  }
  
  pop$transitions$rate <- rate
  # Reculculate total rate.
  pop$totalRate <- sum(pop$transitions$rate)			
  return(pop)
}




#####################################################################################################################

# Calculate initial rates i.e. for all colonies/pathogens. 

#' Calculate birth rates for each colony.
#'
#' Calculate starting birth rate for each colony. Later functions recalculate relevant rates after an event.
#'
#'@inheritParams randEvent
#'@name birthR
#'@family initialRates
#' 
birthR <- function(pop, t){ 
  return(pop$parameters['birth']*colSums(pop$I[,,t]) )
}


#' Calculate death rates for each colony.
#'
#' Calculate starting death rate for each colony. Later functions recalculate relevant rates after an event.
#'
#'@inheritParams randEvent
#'@name deathR
#'@family initialRates
#' 
deathR <- function(pop, t){
  return(as.vector(t(pop$I[, , 1] * pop$parameters['death'] + pop$nInfections * pop$parameters['infectDeath'])))
}

#' Calculate dispersal rates for each colony.
#'
#' Calculate starting dispersal rate for each colony. Later functions recalculate relevant rates after an event.
#'
#'@inheritParams randEvent
#'@name dispersalR
#'@family initialRates
#' 
dispersalR <- function(pop, t){
  return(pop$parameters['dispersal']*as.vector(t(pop$I[, pop$edgeList[,1], t])*pop$adjacency[pop$edgeList]))
}


#' Calculate infection rates for each colony.
#'
#'
#'@inheritParams randEvent
#'@name infectionR
#'@family initialRates
#'@return nColonies * nPathogens vector Grouped by pathogen
#' 

infectionR <- function(pop, t){
  return(as.vector(sapply(1:pop$parameters['nPathogens'], 
    function(rho) pop$parameters['transmission']*pop$I[1, , t] * colSums(pop$I[pop$whichClasses[, rho], , t]))))
}

#' Calculate coinfection rates for each colony.
#'
#'
#'@inheritParams randEvent
#'@name coinfectionR
#'@family initialRates
#'@return nColonies * something. Grouped by to class, then from class
#' 


coinfectionR <- function(pop, t){

  coinfectionTrans <- pop$transitions[pop$transitions$type == 'infection' & pop$transitions$toClass > (pop$parameters['nPathogens'] + 1), ]

  # rate is alpha * beta * fromClass * sumAdditional
  sumAdditions <- sapply(1:length(pop$diseaseAdded), function(i) sum(pop$I[pop$whichClasses[, pop$diseaseAdded[i]], coinfectionTrans$fromColony[i], t]))

  rate <- pop$parameters['transmission'] * pop$parameters['crossImmunity'] * sumAdditions *   pop$I[cbind(coinfectionTrans$fromClass, coinfectionTrans$fromColony, t)]
  return(rate)
}




#' Calculate recovery rates for each colony.
#'
#'
#'@inheritParams randEvent
#'@name infectionR
#'@family initialRates
#'@return nColonies * nPathogens vector Grouped by pathogen
#' 

recoveryR <- function(pop, t){

  # Indices (note reversed order of columns)
  recoveryTrans <- as.matrix(pop$transitions[pop$transitions$type == 'recovery', c(3,2)])
  curPop <- pop$I[, , t]
  return(pop$parameters['recovery'] * curPop[recoveryTrans])
}



#' Find possible possible infection transitions.
#'
#'
#'@inheritParams randEvent
#'@name infectionTrans
#'@family initialRates
#' 
infectionTrans <- function(pop){
  transMatr <- matrix(0, ncol = pop$nClasses, nrow = pop$nClasses)                 
  for(i in 1:pop$nClasses){
    for(j in 1:pop$nClasses){
      setDiff <- length(setdiff(pop$diseaseList[[j]], pop$diseaseList[[i]]))==1
      increasing <- all(pop$diseaseList[[i]] %in% pop$diseaseList[[j]])
      if(setDiff & increasing){transMatr[i,j] <- 1 }
    }
  } 
  return(which(transMatr == 1, arr.ind = TRUE))
}


#' Make vector of which disease is added for each coinfection transitions in pop$transitions
#'
#'
#'@inheritParams randEvent
#'@name findDiseaseAdded
#'@family initialRates
#' 
findDiseaseAdded <- function(pop){

  # What are the possible transitions (col 1 from, col 2 to class.)
  infectTrans <- infectionTrans(pop)


  # Build the infection bit of pop$transisions  
  infectionTransitions <- data.frame(type = 'infection', fromColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])) , fromClass = rep(infectTrans[,1], each = pop$parameters['nColonies']), toColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])), toClass = rep(infectTrans[,2], each = pop$parameters['nColonies']), rate = NA)

  # single infection 'to' class
  to <- infectionTransitions$toClass[infectionTransitions$type == 'infection']
  # just want coinfection 'to' class
  coinfectionTo <- to[to > pop$parameters['nPathogens'] + 1]

  # From classes
  from <- infectionTransitions$fromClass[infectionTransitions$type == 'infection']
  coinfectionFrom <- from[from > 1]

  # which disease is in 
  diseaseAdded <- apply(cbind(coinfectionFrom, coinfectionTo), 1, function(r) pop$diseaseList[[r[2]]][!pop$diseaseList[[r[2]]] %in% pop$diseaseList[[r[1]]]])
  return(diseaseAdded)
}




####################################################################################################################
####################################################################################################################

#' Draw adjacency matrix between colonies
#'
#' For unweighted networks, draw an adjacency matrix showing which colonies are connected and which aren't
#'
#'@inheritParams weightMatrix
#'
#'@name adjMatrix
#'@family Networks 


adjMatrix <- function(locations, maxDistance){
	assert_that(NCOL(locations)==2, NROW(locations)>0, is.numeric(locations))
	assert_that(length(maxDistance)==1, is.numeric(maxDistance))

	d <- as.matrix(dist(locations, diag=TRUE, upper=TRUE))
	a <- (d<maxDistance)*1
	diag(a) <- 0
	
	return(a)
}

#' Draw weighted adjacency matrix between colonies
#'
#' For weighted networks, draw an adjacency matrix showing which colonies are connected and how strongly connected they are
#'   depending on a kernel that defines how connection strength decreases with distance.
#'   
#' Many kernels have a maximum threshold above which colonies aren't connected. This simplified things as otherwise all colonies
#'   would be connected to some extent. This threshold is defined with the argument maxDistance.  
#'
#'@param locations x, y numeric locations of each colony.
#'@param maxDistance For kernels with an upper threshold
#'@param kern Character, giving the name of the distance kernel to be used.
#'
#'@name weightMatrix
#'@family Networks 

weightMatrix <- function(locations, maxDistance, kern){
	assert_that(NCOL(locations)==2, NROW(locations)>0, is.numeric(locations))
	assert_that(length(maxDistance)==1, is.numeric(maxDistance))			

	d <- as.matrix(dist(locations, diag=TRUE, upper=TRUE))
	w <- do.call(kern, list(d, maxDistance))
	return(w)
}

#' Count the number of components in population network
#'
#' If sections of the network are completely unconnected to the rest of the population this has odd effects for an analysis.
#'   Therefore check how many components there are. 
#'   
#'@param W A weighted adjacency matrix
#'
#'@seealso \code{\link{checkIfGiant}}
#'@name numComponents
#'@family Networks 
numComponents <- function(W){
	assert_that(is.matrix(W), NROW(W)==NCOL(W), is.numeric(W))
                    A <- (W > 0)*1
	g <- graph.adjacency(A)
	c <- length(decompose.graph(g))
	return(c)
}

#' Is the network one giant component or multiple small separate networks.
#'
#' If sections of the network are completely unconnected to the rest of the population this has odd effects for an analysis.
#'   Therefore check how many components there are. 
#'   
#'@param A An adjacency matrix
#'
#'@seealso \code{\link{numComponents}}
#'@name checkIfGiant
#'@family Networks 
checkIfGiant <- function(A){
	c <- numComponents(A)
	cig <- if(c==1){ TRUE } else {FALSE}
	return(cig)
}

#' Calculate a dispersal matrix
#'
#' Calulate a weighted (or unweighted) adjacency matrix. Then weight this value to convert to a dispersal matrix.
#'   An element i, j is the proportion of individuals dispersing from i that go to j.
#'   
#'@inheritParams weightMatrix
#'
#'@name popMatrix
#'@family Networks 
popMatrix <- function(locations, maxDistance, kern){
		W <- weightMatrix(locations, maxDistance, paste0(kern, 'Kern')) 
		A <- W/rowSums(W)
	return(A)
}

#' Find the network degree of colonies
#'
#'   
#'@inheritParams seedPathogen
#'
#'@name netDegree
#'@family Networks 
#'@export
netDegree <- function(pop){	
	A <- (pop$adjacency>0)*1
	return(rowSums(A))
}

#' Find the network strength of colonies
#'
#' The strength of a node is the sum of it's weighted edges.
#'   
#'@inheritParams seedPathogen
#'
#'@name netStrength
#'@family Networks 
#'@export
netStrength <- function(pop){
	return(rowSums(pop$adjacency))
}
                


########################################################################################################################

#' Create a list of edges with a column of weights
#'
#'   
#'@inheritParams seedPathogen
#'
#'@name edgeList
#'@family Networks 

edgeList <- function(pop){
	assert_that(is.numeric(pop$locations), NCOL(pop$locations)==2, NROW(pop$locations)>0)

	W <- weightMatrix(pop$locations, pop$parameters['maxDistance'], paste0(pop$models$kernel, 'Kern'))
	listOfEdges <- which(W>0, arr.ind=TRUE)
  
  #add weights
  listOfEdges <- cbind(listOfEdges, W[listOfEdges])
	
	return(listOfEdges)
}





#' Returns an nColony x nPathogens matrix 
#' 
#' To look at global infection, rather than colony-wise infection, this function sums across colonies
#'   to give an nColony x nPathogens for a certain time, t.
#'
#'@inheritParams seedPathogen
#'@param t Time. An integer giving the event number (i.e. the 1000th event).
#'@name sumI

sumI <- function(pop, t){
	return(sapply(1:pop$parameters['nPathogens'], function(x) colSums(pop$I[pop$whichClasses[, x], , t])))
}







