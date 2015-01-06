
#####################################################################################################################

# Calculate rates i.e. for all colonies/pathogens. 


#' Calculate birth rates for each colony.
#'
#' Calculate starting birth rate for each colony. Later functions recalculate relevant rates after an event.
#'
#'@inheritParams randEvent
#'@name birthD
#'@family initialRates
#' 
birthD <- function(pop, t){ 
  return(pop$parameters['birth']*colSums(pop$I[,,t]) )
}


#' Calculate death rates for each colony.
#'
#' Calculate starting death rate for each colony. Later functions recalculate relevant rates after an event.
#'
#'@inheritParams randEvent
#'@name deathD
#'@family initialRates
#' 
deathD <- function(pop, t){
  return(as.vector(t(pop$I[, , t] * (pop$parameters['death'] + pop$nInfections * pop$parameters['infectDeath']))))
}

#' Calculate dispersal rates for each colony.
#'
#' Calculate starting dispersal rate for each colony. Later functions recalculate relevant rates after an event.
#'
#'@inheritParams randEvent
#'@name dispersalD
#'@family initialRates
#' 
dispersalD <- function(pop, t){
  return(pop$parameters['dispersal']*as.vector(t(pop$I[, pop$edgeList[,1], t])*pop$adjacency[pop$edgeList]))
}


#' Calculate infection rates for each colony.
#'
#'
#'@inheritParams randEvent
#'@name infectionD
#'@family initialRates
#'@return nColonies * nPathogens vector Grouped by pathogen
#' 

infectionD <- function(pop, t){
  return(as.vector(sapply(1:pop$parameters['nPathogens'], 
    function(rho) pop$parameters['transmission']*pop$I[1, , t] * colSums(pop$I[pop$whichClasses[, rho], , t]))))
}

#' Calculate coinfection rates for each colony.
#'
#'
#'@inheritParams randEvent
#'@name coinfectionD
#'@family initialRates
#'@return nColonies * something. Grouped by to class, then from class
#' 


coinfectionD <- function(pop, t){

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
#'@name infectionD
#'@family initialRates
#'@return nColonies * nPathogens vector Grouped by pathogen
#' 

recoveryD <- function(pop, t){

  # Indices (note reversed order of columns)
  recoveryTrans <- as.matrix(pop$transitions[pop$transitions$type == 'recovery', c(3,2)])
  curPop <- pop$I[, , t]
  return(pop$parameters['recovery'] * curPop[recoveryTrans])
}




#####################################################################################################################

# Calculate rates for a single colony
# These functions are used for non dispersal events as they are quicker



#' Calculate birth rates for a colony.
#'
#' Calculate starting birth rate for each colony. Later functions recalculate relevant rates after an event.
#'
#'@param pop A population object
#'@param t Integer time step
#'@param colony Integer giving which colony to calculate for.
#'@name birthR
#'@family initialRates
#' 
birthR <- function(pop, t, colony){ 
  return(pop$parameters['birth']*sum(pop$I[, colony, t]) )
}


#' Calculate death rates for a colony.
#'
#' Calculate starting death rate for each colony. Later functions recalculate relevant rates after an event.
#'
#'@inheritParams birthR
#'@name deathR
#'@family initialRates
#' 
deathR <- function(pop, t, colony){
  return(pop$I[, colony, t] * (pop$parameters['death'] + pop$nInfections * pop$parameters['infectDeath']))
}

#' Calculate dispersal rates for a colony.
#'
#' Calculate starting dispersal rate for each colony. Later functions recalculate relevant rates after an event.
#'
#'@inheritParams birthR
#'@name dispersalR
#'@family initialRates
#' 
dispersalR <- function(pop, t, colony){
  return(pop$parameters['dispersal']*as.vector(t(pop$I[, pop$edgeList[pop$edgeList[,2] == colony, 1], t]))*pop$adjacency[pop$edgeList[pop$edgeList[,2] == colony, ]])
}


#' Calculate infection rates for a colony.
#'
#'
#'@inheritParams birthR
#'@name infectionR
#'@family initialRates
#'@return nColonies * nPathogens vector Grouped by pathogen
#' 

infectionR <- function(pop, t, colony){
  as.vector(sapply(1:pop$parameters['nPathogens'], 
    function(rho) pop$parameters['transmission']*pop$I[1, colony, t] * sum(pop$I[pop$whichClasses[, rho], colony, t])))
}

#' Calculate coinfection rates for a colony.
#'
#'
#'@inheritParams birthR
#'@name coinfectionR
#'@family initialRates
#'@return nColonies * something. Grouped by to class, then from class
#' 


coinfectionR <- function(pop, t, colony){

  coinfectionTrans <- pop$transitions[pop$transitions$type == 'infection' & pop$transitions$toClass > (pop$parameters['nPathogens'] + 1) & pop$transitions$fromColony == colony, ]

  oneColDiseaseAdded <- pop$diseaseAdded[seq(1, length(pop$diseaseAdded), pop$parameters['nColonies'])]
  # rate is alpha * beta * fromClass * sumAdditional
  sumAdditions <- sapply(1:length(oneColDiseaseAdded), function(i) sum(pop$I[pop$whichClasses[, oneColDiseaseAdded[i]], coinfectionTrans$fromColony[i], t]))

  rate <- pop$parameters['transmission'] * pop$parameters['crossImmunity'] * sumAdditions *   pop$I[cbind(coinfectionTrans$fromClass, rep(colony, length(coinfectionTrans$fromClass)) , t)]
  return(rate)
}




#' Calculate recovery rates for a colony.
#'
#'
#'@inheritParams birthR
#'@name infectionR
#'@family initialRates
#'@return nColonies * nPathogens vector Grouped by pathogen
#' 

recoveryR <- function(pop, t, colony){

  # Indices (note reversed order of columns)
  recoveryTrans <- as.matrix(pop$transitions[pop$transitions$type == 'recovery' & pop$transitions$fromColony == colony, 3])
  curPop <- pop$I[, colony, t]
  return(pop$parameters['recovery'] * curPop[recoveryTrans])
}

