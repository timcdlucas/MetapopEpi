
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

unweightedKern <- function(distMatrix, thresh){
  W <- (distMatrix<thresh)*1
  diag(W) <- 0
  return(W)
}

