#' Print method for objects of class MetapopEpi
#'
#'@export
#'@name print.MetapopEpi

print.MetapopEpi <- function(p){
  cat('MetapopEpi object.\n\n')
  
  cat('Epidemiological parameters:\n')
  cat('Model: ', p$models$model, '\n')
  print(signif(p$parameters, 2)[c('transmission', 'recovery', 'crossImmunity', 'nPathogens')])
  cat('\n')
  
  cat('Population parameters:\n')
  print(signif(p$parameters, 2)[c('dispersal', 'nColonies', 'maxDistance')])
  cat('\nArray Dimensions:\n', dim(p$I), '\n')
}
  
