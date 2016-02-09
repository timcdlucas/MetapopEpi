
####################################################################################################################
####################################################################################################################



# Population size distributions


#' Create equal colony sizes.
#'
#'@inheritParams exponentialPop 
#'@param colonySize The size of all colonies
#'@name equalPop


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


exponentialPop <- function(nColonies, meanColonySize){
	S <- rexp(nColonies, 1/meanColonySize) %>% ceiling
	return(S)		
}

#' Create colony sizes based on poisson distribution.
#'
#'@inheritParams exponentialPop 
#'@name poissonPop

poissonPop <- function(nColonies, meanColonySize){
	S <- rpois(nColonies, meanColonySize)
	return(S)
}
