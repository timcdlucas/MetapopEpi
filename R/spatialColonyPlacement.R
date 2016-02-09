

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


circleColony <- function(N, space){     
	assert_that(N%%1==0, N>0, is.numeric(space), space>0)
  deg <- 1:N*2*pi/N
  r <- 0.4*space
  x <- cbind(r*cos(deg)+space/2, r*sin(deg)+space/2)
}



