



#'A package containing functions to run multi pathogen, metapopulation epidemiological simulations.
#'
#'@name MetapopEpi
#'@author Tim CD Lucas
#'@docType package
#'@import assertthat RColorBrewer magrittr igraph

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
	is.count(nPathogens)
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
	
	pop$I <- array(0, dim=c(pop$nClasses, nColonies, events), dimnames=c('pathogen', 'colony', 'events'))
	
	# Create colony sizes based on model in colonyDistr
	pop$I[1,,1] <- do.call(paste0(colonyDistr, 'Pop'), list(nColonies, meanColonySize))


	# Make element that stores transition probabilities
	# How many coinfection transitions
	coinfectionTrans <- sum(sapply(0:nPathogens, function(k) choose(nPathogens, k)*(nPathogens-k) ) )

	# How many total transitions
	#distinctTrans <- nColonies + nColonies*pop$nClasses + coinfectionTrans +nColonies*(nColonies-1)*pop$nClasses
	infectTrans <- infectionTrans(pop)
	
	pop$transitions <- rbind(
				data.frame(type='birth', fromColony=NA, fromClass=NA, toColony=1:nColonies, toClass=1, rate=birthR(pop,1)),
				data.frame(type='death', fromColony=rep(1:nColonies, pop$nClasses), fromClass=rep(1:pop$nClasses, each=nColonies), toColony=NA, toClass=NA, rate=deathR(pop,1)),
				data.frame(type='infection', fromColony=rep(1:nColonies, length(infectTrans[,1])) , fromClass=rep(infectTrans[,1], each=nColonies), toColony=rep(1:nColonies, length(infectTrans[,1])), toClass=rep(infectTrans[,2], each=nColonies), rate=NA),
				data.frame(type='dispersal', fromColony=rep(pop$edgeList[,1], pop$nClasses), fromClass=rep(1:pop$nClasses, each=NROW(pop$edgeList)), toColony=rep(pop$edgeList[,2], pop$nClasses), toClass=rep(1:pop$nClasses, each=NROW(pop$edgeList)), rate=NA)
				)

	pop$transitions$rate <- c(birthR(pop,1), deathR(pop,1),  rep(0,nColonies*length(infectTrans[,1])) , dispersalR(pop,1))
                    
  pop$totalRate <- sum(pop$transitions$rate)			

  # Make vector for waiting times.
  pop$waiting <- rep(0,events)
  # Precalculate random uniform numbers. Will go into function waitingTime() to give exponential waiting time
  pop$randU <- runif(events)

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


seedPathogen <- function(pop, pathogens = 1){
	
  assert_that(all(sapply(pathogens, is.count)))
	# can't seed more species than were initiliased.
	assert_that(max(pathogens) <= pop$parameters['nPathogens'])
  assert_that(length(pathogens) <= pop$parameters['nPathogens'])
	
	for(path in pathogens){
		# choose random colony
		r <- sample(pop$parameters['nColonies'], 1)
		pop$I[path+1, r, 1] <- 1
		# keep pop size constant
		pop$I[1, r, 1] <- pop$I[1, r, 1] - 1
	}
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
# @tim Need rewriting. Yup, totally wrong.
waitingTime <- function(randU, S, I, beta, b, d){
	lambda <- I*S*beta + b*(I+S) + d*(I+S)
	w <- log(1-randU)/(-lambda)
	return(w)
}

# 
randEvent <- function(){

        }



# Calc new transition rates.

transRates <- function(pop){

	# Make element that stores transition probabilities
	# How many coinfection transitions
	coinfectionTrans <- sum(sapply(0:pop$parameters['nPathogens'], function(k) choose(pop$parameters['nPathogens'], k)*(pop$parameters['nPathogens']-k) ) )

	# How many total transitions
	#distinctTrans <- nColonies + nColonies*pop$nClasses + coinfectionTrans +nColonies*(nColonies-1)*pop$nClasses
	infectTrans <- infectionTrans(pop)


	transitions <- rbind(
				data.frame(type = 'birth', fromColony = NA, fromClass = NA, toColony = 1:pop$parameters['nColonies'], toClass = 1, rate = birthR(pop,1)),
				data.frame(type = 'death', fromColony = rep(1:pop$parameters['nColonies'], pop$nClasses), fromClass = rep(1:pop$nClasses, each = pop$parameters['nColonies']), toColony = NA, toClass = NA, rate = deathR(pop,1)),
				data.frame(type = 'infection', fromColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])) , fromClass = rep(infectTrans[,1], each = pop$parameters['nColonies']), toColony = rep(1:pop$parameters['nColonies'], length(infectTrans[,1])), toClass = rep(infectTrans[,2], each = pop$parameters['nColonies']), rate = NA),
				data.frame(type = 'dispersal', fromColony = rep(pop$edgeList[,1], pop$nClasses), fromClass = rep(1:pop$nClasses, each = NROW(pop$edgeList)), toColony = rep(pop$edgeList[,2], pop$nClasses), toClass = rep(1:pop$nClasses, each = NROW(pop$edgeList)), rate = NA)
				)

  pop$transitions <- transitions
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
  return(pop$parameters['birth']*(pop$I[1,,t] + colSums(pop$I[,,t])) )
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
  return(as.vector(apply(pop$I[,,t], 1, function(x) x*pop$parameters['death'] )))
}

#' Find possible possible infection transitions.
#'
#'
#'@inheritParams sumI
#'@name infectionTrans
#'@family initialRates
#' 
infectionTrans <- function(pop){
  transMatr <- matrix(0, ncol=pop$nClasses, nrow=pop$nClasses)                 
  for(i in 1:pop$nClasses){
    for(j in 1:pop$nClasses){
      setDiff <- length(setdiff(pop$diseaseList[[j]], pop$diseaseList[[i]]))==1
      increasing <- all(pop$diseaseList[[i]] %in% pop$diseaseList[[j]])
      if(setDiff & increasing){transMatr[i,j] <- 1 }
    }
  } 
  return(which(transMatr==1, arr.ind=TRUE))
}

		
#' Calculate dispersal rates for each colony.
#'
#' Calculate starting dispersal rate for each colony. Later functions recalculate relevant rates after an event.
#'
#'@inheritParams sumI
#'@name dispersalR
#'@family initialRates
#' 
dispersalR <- function(pop,t){
  return(pop$parameters['dispersal']*as.vector(t(pop$I[,pop$edgeList[,1],t])*pop$adjacency[pop$edgeList]))
}

########################################################################################

death <- function(pop, t){
                }

birth <- function(pop, t){
                }

infection <- function(pop, t){
                }

dispersal <- function(pop, t){
                }


####################################################################################################################
####################################################################################################################
# colony placement models




#' Make a matrix of N colony locations placed by a spatial poisson process.
#'
#' Uniformally random placement of colonies in two dimensions.
#'
#'@param N Number of colonies. A single, positive integer.
#'@param space The size of space as given by the length of one side of square space.
#'@name uniformColony
#'@family colony.placement

uniformColony <- function(N, space){
	assert_that(N%%1==0, N>0, is.numeric(space), space>0)
	x <- cbind(runif(N, 0, space), runif(N, 0, space))
	return(x)
}

#' Place colonies in a circle. Most good for viz when all colonies are connected etc..
#'
#' Radial placement of colonies in two dimensions.
#'
#'@inheritParams uniformColony
#'@name circleColony
#'@family colony.placement

circleColony <- function(N, space){     
	assert_that(N%%1==0, N>0, is.numeric(space), space>0)
  deg <- 1:N*2*pi/N
  r <- 0.4*space
  x <- cbind(r*cos(deg)+space/2, r*sin(deg)+space/2)
}






####################################################################################################################
####################################################################################################################



# Population size distributions


#' Create equal colony sizes.
#'
#'@inheritParams exponentialPop 
#'@param colonySize The size of all colonies
#'@name equalPop
#'@family colonysize.distribution

equalPop <- function(nColonies, colonySize){
	S <- rep(colonySize, nColonies)
	return(S)
}

#' Create exponential colony size distribution.
#'
#' Create random colony size distirbution from exponential distribution.
#'
#'@param nColonies The number of colonies. A single integer. 
#'@param meanColonySize Mean colony size. 
#'@name exponentialPop
#'@family colonysize.distribution

exponentialPop <- function(nColonies, meanColonySize){
	S <- rexp(nColonies, 1/meanColonySize) %>% ceiling
	return(S)		
}

#' Create colony sizes based on poisson distribution.
#'
#'@inheritParams exponentialPop 
#'@name poissonPop
#'@family colonysize.distribution
poissonPop <- function(nColonies, meanColonySize){
	S <- rpois(nColonies, meanColonySize)
	return(S)
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
                

####################################################################################################################
####################################################################################################################

#' Inverse kernel for calculating a weight matrix from a distance matrix
#'
#' 1/distance with a threshold. 
#'   
#'@param distMatrix A distance matrix of the euclidean distances between each pair of colonies
#'@param thresh A maximum distance above which colonies are not connected.
#'
#'@name inverseKern
#'@family kernels 

inverseKern <- function(distMatrix, thresh){
	t <- max(distMatrix)/thresh
	W <- max(distMatrix)/distMatrix
	W[is.infinite(W)] <- 0
	W[W<t] <- 0
	return(W)
}

#' Linear kernel for calculating a weight matrix from a distance matrix
#'
#' 1 - distance with a threshold. 
#'   
#'@inheritParams inverseKern
#'
#'@name linearKern
#'@family kernels 

linearKern <- function(distMatrix, thresh){
	t <- max(distMatrix) - thresh
	W <- max(distMatrix) - distMatrix
	W[is.infinite(W)] <- 0
	W[W<t] <- 0
	return(W)
}

#' Unweighted kernel for calculating a weight matrix from a distance matrix
#'
#' Creates an unweighted network. As long as distance < theshold, the weight is 1.
#'   
#'@inheritParams inverseKern
#'
#'@name unweightedKern
#'@family kernels 
unweightedKern <- function(distMatrix, thresh){
  W <- (distMatrix<thresh)*1
  diag(W) <- 0
  return(W)
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







