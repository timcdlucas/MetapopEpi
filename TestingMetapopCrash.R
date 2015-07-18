# 2015-03-24 sims

## The plan
- Examine dispersal with just a few values
- Vary transmission as well
- Quite a few colonies
- Low seed value to mimic evo of new strain
   


```r
library(parallel)
library(MetapopEpi)
library(magrittr)

seed <- 150324

set.seed(seed)

setwd('~/Dropbox/phd/Analysis/pathogenDiversity/Data/150324')
```


```r
# x = 1 to get error at 8%

fullSim <- function(x){

  dispVec <- rep(c(0.001, 0.01, 0.1), 1001/3)
  disp <- dispVec[x]

  tranVec <- rep(c(2, 5, 10), 1001/3)
  tran <- tranVec[x]


  simSeed <- paste0(seed, x)
  set.seed(simSeed)


  #854997 gives error 
  s <- sample(1:1e6, 1)


  # Beginning of runSim

  if (end == 'end'){
    end <- pop$parameters['events']
  }
  
    
  #for (t in start:end){
  # error occurs at t = 40001, but need to run t = 40000 first.
  t = start # = 40000
    

    tMod <- (t - 1) %% pop$parameters['sample'] + 1 #  = 1000


    pop <- randEvent(pop, t, tMod)

    if(tMod == pop$parameters['sample']){
      pop$sampleWaiting[t/pop$parameter s['sample'] + 1] <- sum(pop$waiting)
      pop$sample[, , t/pop$parameters['sample'] + 1] <- pop$I[, , tMod + 1]
      pop$I[, , 1] <- pop$I[, , tMod + 1]
    }
  
# Beginning of randEvent(pop, t, tMod)
 
  event <- pop$transitions[findInterval(pop$randEventU[t] * pop$totalRate, cumsum(c(0, pop$transitions$rate))), ]

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


  s <- 854997
  print(s)
  set.seed(s)
  p <- makePop(model = 'SIR', events = 47000, nColonies = 10, nPathogens = 2,   recovery = 1,  sample = 1000, dispersal = disp, birth = 0.05, death = 0.05, crossImmunity = 0.1, meanColonySize = 3000, infectDeath = 0, transmission = tran, maxDistance = 100)


  for(i in 1:10){
    p <- seedPathogen(p, 2, n = 200, diffCols = FALSE)
  }

  p <- runSim(p, end = 40000)


  pop <- p
  pop$I[2, 1, 40] <- 10
  

  ww <- pop$transitions
  ww$rate <- round(ww$rate, 2)

  pop <- transRates(pop, 1001)

  ww$rate2 <- round(pop$transitions$rate, 2)

  pop <- runSim(pop, start = 40000, end = 'end')







  s <- 854997
  print(s)
  set.seed(s)

  invadeT <- 40000
  sample <- 20

  p <- makePop(model = 'SIR', events = 47000, nColonies = 10, nPathogens = 2,   recovery = 1,  sample = sample, dispersal = disp, birth = 0.05, death = 0.05, crossImmunity = 0.1, meanColonySize = 3000, infectDeath = 0, transmission = tran, maxDistance = 100)


  for(i in 1:10){
    p <- seedPathogen(p, 2, n = 200, diffCols = FALSE)
  }

  p <- runSim(p, end = invadeT)
  p$I[2, 1, (invadeT + 1) %% sample] <- 10
  
  # Must be sample + 1

  p <- transRates(p, (invadeT + 1) %% sample)

  p <- runSim(p, start = invadeT, end = 'end')



















  findDisDistr(p)[1] > 0



  message(paste("finished", x, "Invasion ", invasion ))

  file <- paste0('~/Dropbox/phd/Analysis/pathogenDiversity/Data/150324/150324_', x, '.RData')
  save(p, file = file)
  rm(p)

}
```



d = sapply(2:1001,   function(x) sum(p$I[,,x] - p$I[,,x-1]))



```r
z <- mclapply(101:1000, . %>% fullSim, mc.preschedule = FALSE, mc.cores = 7L)
```

```
## Warning in mclapply(101:1000, . %>% fullSim, mc.preschedule = FALSE,
## mc.cores = 7L): 228 function calls resulted in an error
```

```r
save(z, file = '~/Dropbox/phd/Analysis/pathogenDiversity/Data/150324/errors2.RData')
```

                
