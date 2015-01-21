


#' Make a new population
#'
#' Make a new population at time 0. The resultant population object will then contain all the information needed 
#'   to update the population. 
#'
#'@param model Define the model type with a string e.g. 'SI', 'SIS'.
#'@param nColonies How many colonies/subpopulations do you want.
#'@param colonyDistr The model for colony size distirbution. Currently one of 'equal', 'exponential', 'poisson'. 
#'@param space How large should the space be. A positive number giving the length along a side.
#'@param maxDistance The maximum distance within which two colonies can be connected in the network.
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
#'@param sample The sample rate of population states to store. Default is to store every 1000 events.
#'@param infectDeath Additional death rate due to infection. Death rate will be death + infectDeath * n. infections. 
#'@return A large list that contains the model parameters (numeric in 'parameters' and strings in 'models'. 
#'   The actual time course of the population is in 'I', nColonies x nPathogens x events array. 
#'   The spatial locations of colonies are in 'locations' and 'adjacency' is the weighted (if applicable)
#'   adjacency matric for the population network. 'waiting' is the interevent times. diseaseClasses, diseaseList 
#'   and whichClasses are mostly references for the relationships between multidisease states.
#'@name makePop
#'@family Run.sims
#'@export
#'@examples 
#'p <- makePop()
#'p
#'


makePop <- function(model = 'SIS', nColonies = 5, colonyDistr = 'equal', space = 100, maxDistance = 200, kernel = 'unweighted', events = 10000, colonySpatialDistr = 'uniform', nPathogens = 3, meanColonySize = 10000, birth = 0.001, death = 0.001, dispersal = 0.001, transmission = 0.01, recovery = 0.005, crossImmunity = 0.1, sample = 1000, infectDeath = 0){
			
	# Define lists of available options
	modelList <-  c('SI', 'SIS', 'DISEASEFREE', 'SIR')
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
  assert_that(events > sample)

  assert_that(all(is.numeric(c(birth, death, dispersal, transmission, recovery))))
  assert_that(all(c(birth, death, dispersal, transmission, recovery) >= 0))

  assert_that(meanColonySize > 0)

  assert_that(crossImmunity >= 0, crossImmunity <= 1)

	# Now make the population object
	pop <- list(parameters = c(nColonies = nColonies, 
						    space = space, 
						    maxDistance = maxDistance, 
						    events = events, 
						    nPathogens = nPathogens, 
						    meanColonySize = meanColonySize,
                birth = birth,
                death = death,
                transmission = transmission,
                recovery = recovery,
                dispersal = dispersal, 
                crossImmunity = crossImmunity,
                sample = sample, 
                infectDeath = infectDeath),
			    models = data.frame(model = model, 
						    colonyDistr = colonyDistr, 
						    kernel = kernel,
                colonySpatialDistr = colonySpatialDistr,  
                stringsAsFactors = FALSE),
			    locations = NULL,
			    adjacency = NULL,
			    I = NULL)
	
	# Place colonies based on model in colonySpatialDistr
	pop$locations <- do.call(paste0(colonySpatialDistr, 'Colony'), list(nColonies, space))



	# Create adjacency matrix
	pop$adjacency <- popMatrix(pop$locations, maxDistance, kernel)
	pop$edgeList <- which(pop$adjacency != 0, arr.ind = TRUE)
	# Check that the network is connected
	if (!checkIfGiant(pop$adjacency)) warning('The network is made up of multiple unconnected components')			


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
			strings[i] <- paste0(as.character(set[,i]), collapse = '')
			lists[[i]] <- set[, i]
		}
		pop$diseaseClasses <- c(pop$diseaseClasses, strings)
		pop$diseaseList <- c(pop$diseaseList, lists)
	}
	
	# A matrix of which classes contain inds infected with each disease.
	# So all infectives of disease 2 would be sum(pop$I[whichClasses[,2], , t]) or something.
	pop$whichClasses <- sapply(1:nPathogens, function(x) grep(x, pop$diseaseClasses))
	pop$nClasses <- length(pop$diseaseClasses)
  if(pop$models$model == 'SIR') pop$nClasses <- length(pop$diseaseClasses) + 1
  pop$nInfections <- sapply(pop$diseaseList, length)
  if(pop$models$model == 'SIR') pop$nInfections <- c(pop$nInfections, 0)

  # A list of which disease is added in order to reach new class
  #   for all coinfection rows in pop$transitions
  pop$diseaseAdded <- findDiseaseAdded(pop)
	
	pop$I <- array(0, dim = c(pop$nClasses, nColonies, sample + 1), dimnames = c('pathogen', 'colony', 'events'))
	
	# Create colony sizes based on model in colonyDistr
	pop$I[1, , 1] <- do.call(paste0(colonyDistr, 'Pop'), list(nColonies, meanColonySize))


  # Make tmp array.
  pop$sample <- array(0, dim = c(pop$nClasses, nColonies, events/sample + 1), dimnames = c('pathogen', 'colony', 'events'))
	pop$sample[1,,1] <- do.call(paste0(colonyDistr, 'Pop'), list(nColonies, meanColonySize))

  pop <- initTransitions(pop)
  # Calculate transition rates. Set dispersal = TRUE so that ALL rates are calculated, not just in a given colony
	pop <- transRates(pop, 1, 'dispersal')

  # Make vector for waiting times both sample and final
  pop$waiting <- rep(0, sample + 1)
  pop$sampleWaiting <- rep(0, events/sample + 1)
  
  # Precalculate random uniform numbers. Will go into function waitingTime() to give exponential waiting time
  pop$randU <- runif(events + 1)

  # Precalculate random uniform numbers. Will go into function randomEvent() to select an event
  pop$randEventU <- runif(events + 1)

  class(pop) <- 'MetapopEpi'

	return(pop)
}










#' Build the transition data.frame which will have rates populated later.
#'
#'
#'@inheritParams randEvent
#'@name initTransitions
#'@family initialRates
#'@export 



initTransitions <- function(pop){
 
  coinfectionTrans <- sum(sapply(0:pop$parameters['nPathogens'], 
    function(k) choose(pop$parameters['nPathogens'], k) * (pop$parameters['nPathogens'] - k)))


  infectTrans <- infectionTrans(pop)

  pop$transitions <- rbind(
			  data.frame(type = 'birth', fromColony = NA, fromClass = NA, toColony = 1:pop$parameters['nColonies'], toClass = 1, rate = birthD(pop,1), stringsAsFactors = FALSE),
			  data.frame(type = 'death', fromColony = rep(1:pop$parameters['nColonies'], pop$nClasses), fromClass = rep(1:pop$nClasses, each = pop$parameters['nColonies']), toColony = NA, toClass = NA, rate = deathD(pop,1)),
			  data.frame(type = 'infection', fromColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])) , fromClass = rep(infectTrans[,1], each = pop$parameters['nColonies']), toColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])), toClass = rep(infectTrans[,2], each = pop$parameters['nColonies']), rate = NA),
			  data.frame(type = 'dispersal', fromColony = rep(pop$edgeList[,1], pop$nClasses), fromClass = rep(1:pop$nClasses, each = NROW(pop$edgeList)), toColony = rep(pop$edgeList[,2], pop$nClasses), toClass = rep(1:pop$nClasses, each = NROW(pop$edgeList)), rate = NA)
			  )


  if(pop$models$model == 'SIS'){
    # infection in reverse
    pop$transitions <- rbind(pop$transitions,
      data.frame(type = 'recovery', 
                 fromColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])) , 
                 fromClass = rep(infectTrans[, 2], each = pop$parameters['nColonies']), 
                 toColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])), 
                 toClass = rep(infectTrans[, 1], each = pop$parameters['nColonies']), 
                 rate = NA)
    )
  }


  if(pop$models$model == 'SIR'){
    # infection in reverse
    pop$transitions <- rbind(pop$transitions,
      data.frame(type = 'recovery', 
                 fromColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])) , 
                 fromClass = rep(infectTrans[, 2], each = pop$parameters['nColonies']), 
                 toColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])), 
                 toClass = pop$nClasses, 
                 rate = NA)
    )
  }
  return(pop)
}
  
  

