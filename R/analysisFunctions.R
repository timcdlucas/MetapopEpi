

#' Calculate mean disease distributions (how many individuals infected with each disease) at the end of a simulation
#'
#'@param pop A population object from makePop()
#'@param final Integer giving the number of recorded population states to average over. Note that this is the number of recorded events, not actual events so will cover a period of final * sample.
#'@name findDisDistr
#'@export

findDisDistr <- function(pop, final = 10){
  is.count(final)
  end <- pop$parameters['events']/pop$parameters['sample'] + 1
  finalI <- apply(pop$sample[, , (end - final):end], 1, mean)

  finalDis <- sapply(1:pop$parameters['nPathogens'], function(x) finalI[pop$whichClasses[, x]] %>% sum)
}

#' Calculate coinfection disease distributions (how many individuals infected with number of diseases) at the end of a simulation
#'
#'@inheritParams findDisDistr
#'@name findCoinfDistr
#'@export

findCoinfDistr <- function(pop, final = 10){
  is.count(final)
  end <- pop$parameters['events']/pop$parameters['sample'] + 1
  finalI <- apply(pop$sample[, , (end - final):end], 1, mean)

  coinf <- sapply(pop$diseaseList, length)
  
  finalCoinf <- aggregate(finalI, by = list(coinf), FUN = sum)[, 2]

}

