



#'A package containing functions to run multi pathogen, metapopulation epidemiological simulations.
#'
#'@name MetapopEpi
#'@author Tim CD Lucas
#'@docType package
#'@import assertthat RColorBrewer magrittr igraph ggplot2 reshape2 compiler

NULL





#' Make a new population
#'
#' Make a new population at time 0. The resultant population object will then contain all the information needed 
#'   to update the population. 
#'
#'@param model Define the model type with a string e.g. 'SI', 'SIS'.
#'@param nColonies How many colonies/subpopulations do you want.
#'@param colonyDistr The model for colony size distirbution. Currently one of 'equal', 'exponential', 'poisson'. 
#'@param space How large should the space be. A positive number giving the length along a side.
#'@param maxDistace The maximum distance within which two colonies can be connected in the network.
#'@param kernel How should the network be weighted with respect to distance between colonies. 'inverse', 'linear', 'unweighted'
#'@param events Integer defining the number of events to run before stopping the simulation. A single event is an infection, death, birth, migration etc.
#'@param colonySpatialDistr How are the colonies distributed in space. 'uniform', 'circle'
#'@param nPathogens Integer defining the number of pathogen species or strains.
#'@param meanColonySize Colony size.
#'@param birth Positive numeric. Birth rate per individual per unit time.
#'@param death Positive numeric. Death rate per individual per unit time.
#'@param dispersal Positive numeric. Dispersal rate per individual per unit time.
#'@param transmission Positive numeric. Transmission rate per unit time.
#'@param recovery Positive numeric. Recovery rate per individual per unit time. Disease duration is therefore 1/recovery
#'@param crossImmunity Numeric between 0 and 1 that governs cross immunity. 0 is full cross immunity, 1 is no cross immunity.
#'@return A large list that contains the model parameters (numeric in 'parameters' and strings in 'models'. 
#'   The actual time course of the population is in 'I', nColonies x nPathogens x events array. 
#'   The spatial locations of colonies are in 'locations' and 'adjacency' is the weighted (if applicable)
#'   adjacency matric for the population network. 'waiting' is the interevent times. diseaseClasses, diseaseList 
#'   and whichClasses are mostly references for the relationships between multidisease states.
#'@name makePop
#'@family Run.sims
#'@export
#'@examples 
#'p <- buildPop()
#'p
#'


makePop <- function(model = 'SI', nColonies = 5, colonyDistr = 'equal', space = 100, maxDistance = 200, kernel = 'unweighted', events = 1000, colonySpatialDistr = 'uniform', nPathogens = 3, meanColonySize = 10000, birth = 0.001, death = 0.001, dispersal = 0.001, transmission = 0.01, recovery = 0.005, crossImmunity = 0.1){
			
	# Define lists of available options
	modelList <-  c('SI', 'SIS', 'DISEASEFREE')
	colonyDistrList <- c('equal', 'exponential', 'poisson')
	kernelList <- c('inverse', 'linear', 'unweighted')
	colonySpatialDistrList <- c('uniform', 'circle')

	# Test if options are valid
	model <- toupper(model)
	if(!model %in% modelList)
		stop(paste(model,"is not a currently implemented model. Try one of ", paste(modelList, collapse=', ')))
	colonyDistr <- tolower(colonyDistr)
	if(!colonyDistr %in% colonyDistrList)
		stop(paste("Argument colonyDistr='", colonyDistr,"' is not a currently implemented colony size distribution. Try one of ", paste(colonyDistrList, collapse=', ')))
	kernel <- tolower(kernel)
	if(!kernel %in% kernelList)
		stop(paste(kernel,"is not a currently implemented spatial kernel. Try one of ", paste(kernelList, collapse=', ')))
	colonySpatialDistr <- tolower(colonySpatialDistr)
	if(!colonySpatialDistr %in% colonySpatialDistrList)
		stop(paste(colonySpatialDistr,"is not an implemented spatial colony distribution. Try one of ", paste(colonySpatialDistrList, collapse=', ')))
	
	# Some more assertions
	is.count(nColonies)
	is.count(events)
  assert_that(events > 0)
	is.count(nPathogens)
  assert_that(nPathogens > 1)
  assert_that(nColonies > 1)
	assert_that(is.numeric(space), length(space)==1)
	assert_that(is.numeric(maxDistance), length(maxDistance)==1)

  assert_that(all(is.numeric(c(birth, death, dispersal, transmission, recovery))))
  assert_that(all(c(birth, death, dispersal, transmission, recovery) >= 0))

  assert_that(meanColonySize > 0)

  assert_that(crossImmunity >= 0, crossImmunity <= 1)

	# Now make the population object
	pop <- list(parameters=c(nColonies=nColonies, 
						    space=space, 
						    maxDistance=maxDistance, 
						    events=events, 
						    nPathogens=nPathogens, 
						    meanColonySize=meanColonySize,
                birth=birth,
                death=death,
                transmission=transmission,
                recovery=recovery,
                dispersal=dispersal, 
                crossImmunity=crossImmunity),
			    models=data.frame(model=model, 
						    colonyDistr=colonyDistr, 
						    kernel=kernel,
                colonySpatialDistr=colonySpatialDistr,                            
                stringsAsFactors=FALSE),
			    locations=NULL,
			    adjacency=NULL,
			    I=NULL)
	
	# Place colonies based on model in colonySpatialDistr
	pop$locations <- do.call(paste0(colonySpatialDistr, 'Colony'), list(nColonies, space))



	# Create adjacency matrix
	pop$adjacency <- popMatrix(pop$locations, maxDistance, kernel)
	pop$edgeList <- which(pop$adjacency != 0, arr.ind=TRUE)
	# Check that the network is connected
	if(!checkIfGiant(pop$adjacency)) warning('The network is made up of multiple unconnected components')			


	# Create list elements for other disease classed
	# This is for total cross immunity. 
	# if(model %in% c('SIR')) pop$R <- matrix(0,nColonies, events) 
	
	# Now make a list element for each pathogen


	pop$diseaseClasses <- NULL
	pop$diseaseList <- NULL
	
			
	for(k in 0:nPathogens){
		set <- combn(1:nPathogens, k)
		strings <- rep('', NCOL(set))
		lists <- replicate(NCOL(set), list())
		for(i in 1:NCOL(set)){
			strings[i] <- paste0(as.character(set[,i]), collapse='')
			lists[[i]] <- set[,i]
		}
		pop$diseaseClasses <- c(pop$diseaseClasses, strings)
		pop$diseaseList <- c(pop$diseaseList, lists)
	}
	
	# A matrix of which classes contain inds infected with each disease.
	# So all infectives of disease 2 would be sum(pop$I[whichClasses[,2], , t]) or something.
	pop$whichClasses <- sapply(1:nPathogens, function(x) grep(x, pop$diseaseClasses))
	pop$nClasses <- length(pop$diseaseClasses)

  # A list of which disease is added in order to reach new class
  #   for all coinfection rows in pop$transitions
  pop$diseaseAdded <- findDiseaseAdded(pop)
	
	pop$I <- array(0, dim = c(pop$nClasses, nColonies, events + 1), dimnames = c('pathogen', 'colony', 'events'))
	
	# Create colony sizes based on model in colonyDistr
	pop$I[1,,1] <- do.call(paste0(colonyDistr, 'Pop'), list(nColonies, meanColonySize))

  pop <- initTransitions(pop)
	pop <- transRates(pop, 1)

  # Make vector for waiting times.
  pop$waiting <- rep(0,events + 1)
  
  # Precalculate random uniform numbers. Will go into function waitingTime() to give exponential waiting time
  pop$randU <- runif(events + 1)

  # Precalculate random uniform numbers. Will go into function randomEvent() to select an event
  pop$randEventU <- runif(events + 1)

  class(pop) <- 'MetapopEpi'

	return(pop)
}


#' Add single individuals infected with a number of pathogens
#' 
#' Allows addition of an infected individual of 1 or more pathogens. 
#'
#'@param pop A metapopulation object as created by makePop().
#'@param pathogens Character vector of which pathogens should be seeded
#'
#'@return An update object of same structure as from makePop
#'@family Run,sims
#'
#'@export
#'@name seedPathogen
#'@examples
#'p <- makePop()
#'u <- seedPathogen(p, pathogens = 1:p$parameters['nPathogens'])


seedPathogen <- function(pop, pathogens = 1, n = 1){
	
  assert_that(all(sapply(pathogens, is.count)))
	# can't seed more species than were initiliased.
	assert_that(max(pathogens) <= pop$parameters['nPathogens'])
  assert_that(length(pathogens) <= pop$parameters['nPathogens'])
	
	for(path in pathogens){
		# choose random colony
		r <- sample(pop$parameters['nColonies'], 1)
		pop$I[path+1, r, 1] <- n
		# keep pop size constant
		pop$I[1, r, 1] <- pop$I[1, r, 1] - n
	}
  pop <- transRates(pop, 1)

	return(pop)
}


#' Perform a random event
#' 
#' Randomly selects an event weighted by their rates, runs the event and then calculates a waiting time.
#'@param pop A population object
#'@param t Time step. Note the population should be AT time t, going to t+1.
#'@name randEvent 

randEvent <- function(pop, t){
  # Calculate random number each time. Expect, precalcing randUnifs might be quicker.
  # event <- sample(1:nrow(pop$transitions), size = 1, prob = pop$transitions$rate)

  # Using precalcuated random uniforms in [0,1], pop$randEventU
  # Find which interval this falls in as quick way to select value
  # Should have a think about open and closed intervalsfromClass toColony toClass rate

  event <- pop$transitions[findInterval(pop$randEventU[t]*pop$totalRate, cumsum(c(0, pop$transitions$rate))), ]

  # Copy pop to t + 1
  pop$I[, , t + 1] <- pop$I[, , t]
  
  # run event
  pop$I[event[['fromClass']], event[['fromColony']], t + 1] <- pop$I[event[['fromClass']], event[['fromColony']], t ] - 1
  pop$I[event[['toClass']], event[['toColony']], t + 1] <- pop$I[event[['toClass']], event[['toColony']], t] + 1

  pop <- transRates(pop, t + 1)
  
  pop <- waitingTime(pop, t)

  return(pop)
}


  

#' Run simulation up to certain event number
#'
#' Each time step, pick a random event and perform that event.

#'@export
#'@name runSim
#'@inheritParams randEvent
#'@param t The event number to run the simulation until.

runSim <- function(pop, time = 'end'){

  assert_that(is.count(time) | time == 'end')

  if (time == 'end'){
    time <- pop$parameters['events']
  }
  

  pb <- txtProgressBar(1, pop$parameters['events'], style = 3)
    
  for (t in 1:pop$parameters['events']){
    pop <- randEvent(pop, t)
    setTxtProgressBar(pb, t)
  }

  close(pb)
  return(pop)
}

###############################################################################

#' Returns an nColony x nPathogens matrix 
#' 
#' To look at global infection, rather than colony-wise infection, this function sums across colonies
#'   to give an nColony x nPathogens for a certain time, t.
#'
#'@inheritParams seedPathogen
#'@param t Time. An integer giving the event number (i.e. the 1000th event).
#'@name sumI
sumI <- function(pop, t){
	return(sapply(1:pop$parameters['nPathogens'], function(x) colSums(pop$I[pop$whichClasses[,x],,t])))
}




#' Takes a number [0,1] (which is presumably random) and returns waiting time.
#'
#' A single draw from a exponential distribution given randU, a number between 0 and 1.
#'   As population now has 
#'
#'@param randU Numeric in [0,1]. Typically this will a random uniform number. 
#'@name waitingTime

waitingTime <- function(pop, t){
  pop$waiting[t + 1] <- log(1 - pop$randU[t])/(-pop$totalRate)
  return(pop)
}






#' Calculate new transition rates.
#'
#'
#'@inheritParams sumI
#'@name transRates
#'@family initialRates
#'@export 


transRates <- function(pop, t){

  # Infection from which class to which class.

  transitions$rate <- c(birthR(pop, t), deathR(pop, t),  infectionR(pop, t), coinfectionR(pop, t), dispersalR(pop, t))

  pop$transitions <- transitions

  # Reculculate total rate.
  pop$totalRate <- sum(pop$transitions$rate)			
  return(pop)
}




#' Build the transition data.frame which will have rates populated later.
#'
#'
#'@inheritParams sumI
#'@name initTransitions
#'@family initialRates
#'@export 



initTransitions <- function(pop){
 
  coinfectionTrans <- sum(sapply(0:pop$parameters['nPathogens'], function(k) choose(pop$parameters['nPathogens'], k)*(pop$parameters['nPathogens']-k) ) )
  infectTrans <- infectionTrans(pop)

  pop$transitions <- rbind(
			  data.frame(type = 'birth', fromColony = NA, fromClass = NA, toColony = 1:pop$parameters['nColonies'], toClass = 1, rate = birthR(pop,1), stringsAsFactors = FALSE),
			  data.frame(type = 'death', fromColony = rep(1:pop$parameters['nColonies'], pop$nClasses), fromClass = rep(1:pop$nClasses, each = pop$parameters['nColonies']), toColony = NA, toClass = NA, rate = deathR(pop,1)),
			  data.frame(type = 'infection', fromColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])) , fromClass = rep(infectTrans[,1], each = pop$parameters['nColonies']), toColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])), toClass = rep(infectTrans[,2], each = pop$parameters['nColonies']), rate = NA),
			  data.frame(type = 'dispersal', fromColony = rep(pop$edgeList[,1], pop$nClasses), fromClass = rep(1:pop$nClasses, each = NROW(pop$edgeList)), toColony = rep(pop$edgeList[,2], pop$nClasses), toClass = rep(1:pop$nClasses, each = NROW(pop$edgeList)), rate = NA)
			  )
  return(pop)
}
  
  

#####################################################################################################################

# Calculate initial rates i.e. for all colonies/pathogens. 

#' Calculate birth rates for each colony.
#'
#' Calculate starting birth rate for each colony. Later functions recalculate relevant rates after an event.
#'
#'@inheritParams sumI
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
#'@inheritParams sumI
#'@name deathR
#'@family initialRates
#' 
deathR <- function(pop, t){
  return(as.vector(apply(pop$I[,, t], 1, function(x) x*pop$parameters['death'] )))
}

#' Calculate dispersal rates for each colony.
#'
#' Calculate starting dispersal rate for each colony. Later functions recalculate relevant rates after an event.
#'
#'@inheritParams sumI
#'@name dispersalR
#'@family initialRates
#' 
dispersalR <- function(pop, t){
  return(pop$parameters['dispersal']*as.vector(t(pop$I[, pop$edgeList[,1], t])*pop$adjacency[pop$edgeList]))
}


#' Calculate infection rates for each colony.
#'
#'
#'@inheritParams sumI
#'@name infectionR
#'@family initialRates
#'@return nColonies * nPathogens vector Grouped by pathogen
#' 

infectionR <- function(pop, t, infectTrans){
  return(as.vector(sapply(1:pop$parameters['nPathogens'], 
    function(rho) pop$parameters['transmission']*pop$I[1, , t] * colSums(pop$I[pop$whichClasses[, rho], , t]))))
}

#' Calculate coinfection rates for each colony.
#'
#'
#'@inheritParams sumI
#'@name coinfectionR
#'@family initialRates
#'@return nColonies * something. Grouped by to class, then from class
#' 


coinfectionR <- function(pop, t){

  coinfectionTrans <- pop$transitions[pop$transitions$type == 'infection' & pop$transitions$toClass > (pop$parameters['nPathogens'] + 1), ]

  # rate is alpha * beta * fromClass * sumAdditional
  sumAdditions <- colSums(pop$I[pop$whichClasses[, pop$diseaseAdded[1:length(pop$diseaseAdded)]], coinfectionTrans$fromColony[1:length(pop$diseaseAdded)], t])

  rate <- pop$parameters['transmission'] * pop$parameters['crossImmunity'] * sumAdditions *   pop$I[cbind(coinfectionTrans$fromClass, coinfectionTrans$fromColony, t)]
}



#' Find possible possible infection transitions.
#'
#'
#'@inheritParams sumI
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
#'@inheritParams sumI
#'@name findDiseaseAdded
#'@family initialRates
#' 
findDiseaseAdded <- function(pop){
  infectTrans <- infectionTrans(pop)
  toClass <- infectTrans[(pop$parameters['nPathogens'] + 1):(pop$nClasses + pop$parameters['nPathogens'] + 1), 2]


  infectionTransitions <- data.frame(type = 'infection', fromColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])) , fromClass = rep(infectTrans[,1], each = pop$parameters['nColonies']), toColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])), toClass = rep(infectTrans[,2], each = pop$parameters['nColonies']), rate = NA)

  to <- infectionTransitions$toClass[infectionTransitions$type == 'infection']
  coinfectionTo <- to[to > pop$parameters['nPathogens'] + 1]

  from <- infectionTransitions$fromClass[infectionTransitions$type == 'infection']
  coinfectionFrom <- from[from > 1]

  diseaseAdded <- apply(cbind(coinfectionFrom, coinfectionTo), 1, function(r) pop$diseaseList[[r[2]]][!pop$diseaseList[[r[2]]] %in% pop$diseaseList[[r[1]]]])
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







