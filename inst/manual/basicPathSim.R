library(MetapopEpi)


system.time({
pop <- makePop(events = 2000)

pop <- seedPathogen(pop, 1)


pop <- runSim(pop)
}
)



save(pop, file = '~/Dropbox/phd/Analysis/pathogenDiversity/20000events.RData')


load('~/Dropbox/phd/Analysis/pathogenDiversity/Data/20000events.RData')






pop2 <- makePop(events = 200000)

pop2 <- seedPathogen(pop2, 1)



for (t in 1:pop$parameters['events']){
  pop2 <- MetapopEpi:::randEvent(pop2, t)
}
save(pop2)

save(pop2, file = '~/Dropbox/phd/Analysis/pathogenDiversity/200000events.RData')



pAll(pop)

pSus(pop)


aw

d <- makePop()
d <- seedPathogen(d, c(1:3))

Rprof('~/Dropbox/phd/Analysis/pathogenDiversity/Data/prof.out', line.profiling = TRUE)
d2 <- runSim(d)
Rprof(NULL)


summaryRprof('~/Dropbox/phd/Analysis/pathogenDiversity/Data/prof.out')


Rprof('~/Dropbox/phd/Analysis/pathogenDiversity/Data/transTratesprof.out')
for(x in 1:1000) z <- transRates(pop,1)
Rprof(NULL)


summaryRprof('~/Dropbox/phd/Analysis/pathogenDiversity/Data/transTratesprof.out')






pD <- makePop(transmission = 1)
pD <- seedPathogen(pD, 1,  n=1000)
pD <- runSim(pD)
pAll(pD)







t1 <- system.time({
pD <- makePop(events = 1000)
pD <- seedPathogen(pD, c(1:3))
pD <- runSim(pD)
})

pE <- makePop(events = 500000)
pE <- seedPathogen(pE, c(1:3))
pE <- runSim(pE)



m <- microbenchmark(
  coinfection = pop$transitions[pop$transitions$type == 'infection' & pop$transitions$toClass > (pop$parameters['nPathogens'] + 1), ],
  sumAdditions = sapply(1:length(pop$diseaseAdded), function(i) sum(pop$I[pop$whichClasses[, pop$diseaseAdded[i]], coinfectionTrans$fromColony[i], t])),
  rate = pop$parameters['transmission'] * pop$parameters['crossImmunity'] * sumAdditions *   pop$I[cbind(coinfectionTrans$fromClass, coinfectionTrans$fromColony, t)],
  times = 500L
)
autoplot(m, log = T)


