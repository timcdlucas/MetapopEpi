context('Test internal functions for makePop and runnings sims.')


test_that('randEvent works correctly', {
 
  set.seed(2)
  pop <- makePop(model = 'SI', nColonies = 2)
  pop <- seedPathogen(pop, c(1:3), diffCols = FALSE)
  
  # Artificially force t 1:4 to be birth, death, infection, dispersal
  pop$randEventU[1] <- 0.0000001
  pop$randEventU[2] <- 0.1
  pop$randEventU[3] <- 0.5
  pop$randEventU[4] <- 0.999999


  # birth => Susc row one bigger
  pop <- randEvent(pop, 1)

  # One has increased by one, rest haven't changed.
  expect_true(sum(pop$I[1, , 2] - pop$I[1, , 1]) == 1)
  expect_true(sum(pop$I[1, , 2] - pop$I[1, , 1] != 0) == 1)


  # death => whole time slice less one
  pop <- randEvent(pop, 2)

  expect_true(sum(pop$I[, , 3] - pop$I[, , 2]) == -1)
  expect_true(sum(pop$I[1, , 3] - pop$I[1, , 2] != 0) == 1)

  # infection
  pop <- randEvent(pop, 3)
  expect_true(sum(pop$I[, , 4] - pop$I[, , 3]) == 0)
  expect_true(sum(pop$I[, , 4] - pop$I[, , 3] != 0) == 2)
  expect_true(sum(pop$I[, , 4] - pop$I[, , 3] == 1) == 1)
  expect_true(sum(pop$I[, , 4] - pop$I[, , 3] == -1) == 1)


  newi <- which(pop$I[, , 4] - pop$I[, , 3] == 1, arr.ind = TRUE)
  oldi <- which(pop$I[, , 4] - pop$I[, , 3] == -1, arr.ind = TRUE)

  # infection should be in same colony
  expect_true(newi[2] == oldi[2])

  # new infectious class should have more disease than old
  #   Would be better to do pop$diseaseList[[oldi[1]]] %in% pop$diseaseList[[newi[1]]] 
  #   but can't work out how to deal with empty class of Susceptible class.
  expect_true(newi[1] > oldi[1])

  # dispersal one column +1, other -1
  pop <- randEvent(pop, 4)
  expect_true(sum(pop$I[, , 5] - pop$I[, , 4]) == 0)
  expect_true(sum(pop$I[, , 5] - pop$I[, , 4] != 0) == 2)
  expect_true(sum(pop$I[, , 5] - pop$I[, , 4] == 1) == 1)
  expect_true(sum(pop$I[, , 5] - pop$I[, , 4] == -1) == 1)

  newd <- which(pop$I[, , 5] - pop$I[, , 4] == 1, arr.ind = TRUE)
  oldd <- which(pop$I[, , 5] - pop$I[, , 4] == -1, arr.ind = TRUE)

  # dispersal should be between colonies
  expect_true(newd[2] != oldd[2])
  
  # infection class should remain unchanged
  expect_true(newd[1] == oldd[1])

})


test_that('waitingTime works correctly', {
  pop <- makePop(model = 'SI', events = 20, sample = 1)
  
  p1 <- runSim(pop)

  # All waiting times should be positive
  expect_true(all(p1$waiting[-1] > 0))


  # Check that the only thing changed is pop$waiting[t +1]
  pop2 <- makePop(model = 'SI', events = 20, sample = 1)

  pop3 <- waitingTime(pop2, 1)

  restOfPop2 <- pop2
  restOfPop2$waiting <- NULL

  restOfPop3 <- pop3
  restOfPop3$waiting <- NULL

  # Everything except pop$waiting
  expect_equal(restOfPop2, restOfPop3)

  # pop$waiting except t + 1 
  expect_true(all(pop3$waiting[-2] == 0))

})



test_that('transRates works', {
  pop <- makePop(model = 'SI', events = 20, sample = 1)

  pop2 <- transRates(pop, 1)

  expect_equal(pop$transitions, pop2$transitions)

  expect_true(all(pop2$transitions$rate[pop2$transitions$type == 'birth'] == 10))
  expect_true(all(pop2$transitions$rate[pop2$transitions$type == 'death' & pop2$transitions$fromClass == 1] == 10))
  expect_true(all(pop2$transitions$rate[pop2$transitions$type == 'death' & pop2$transitions$fromClass != 1] == 0))
  expect_true(all(pop2$transitions$rate[pop2$transitions$type == 'infection'] == 0))
  expect_true(all(pop2$transitions$rate[pop2$transitions$type == 'dispersal' & pop2$transitions$fromClass == 1] == 2.5))
  expect_true(all(pop2$transitions$rate[pop2$transitions$type == 'dispersal' & pop2$transitions$fromClass != 1] == 0))


  pop3 <- makePop(model = 'SI', nColonies = 2, nPathogens = 3, transmission = 1, dispersal = 2, birth = 3, death = 4)

  pop3$I[, , 1] <- c(1:16)

  pop4 <- transRates(pop3, 1)

  # Dispersal, population * dispersal rate
  expect_true(all(pop4$transitions$rate[pop4$transitions$type == 'dispersal' & pop4$transitions$fromColony == 1] == 1:8 * 2))
  expect_true(all(pop4$transitions$rate[pop4$transitions$type == 'dispersal' & pop4$transitions$fromColony == 2] == 9:16 * 2))

  # birth, sum of colony * birth

  expect_true(pop4$transitions[1, 'rate'] == sum(1:8) * 3)
  expect_true(pop4$transitions[2, 'rate'] == sum(9:16) * 3)

  # death, class * death rate

  expect_true(all(pop4$transitions$rate[pop4$transitions$type == 'death' & pop4$transitions$fromColony == 1] == 1:8 * 4))
  expect_true(all(pop4$transitions$rate[pop4$transitions$type == 'death' & pop4$transitions$fromColony == 2] == 9:16 * 4))
  

  pop5 <- makePop(model = 'SI', nColonies = 2, nPathogens = 3, transmission = 1, crossImmunity = 1, meanColonySize = 1)

  pop6 <- pop5
  pop6$I[1:2, 1, 1] <- 1

  pop6Trans <- transRates(pop6, 1)

  inf <- pop6Trans$transitions$rate[pop6Trans$transitions$fromClass == 1 & pop6Trans$transitions$toClass == 2 & pop6Trans$transitions$toColony == 1 ] 
  inf <- inf[!is.na(inf)]


  expect_equal( inf, 1 * 1 * 1 )


  # Check other classes get added properly.
  pop7 <- pop5
  pop7$I[c(1:2, 5), 1, 1] <- 1

  pop7Trans <- transRates(pop7, 1)

  inf <- pop7Trans$transitions$rate[pop7Trans$transitions$fromClass == 1 & pop7Trans$transitions$toClass == 2 & pop7Trans$transitions$toColony == 1 ] 
  inf <- inf[!is.na(inf)]


  expect_equal( inf, 2 * 1 * 1 )

  
  # coinfection
  pop8 <- pop5
  pop8$I[c(1:3), 1, 1] <- c(0, 1, 1)

  pop8Trans <- transRates(pop8, 1)

  # 1 -> 1 and 2 (i.e. class 2 to class 5)
  inf <- pop8Trans$transitions$rate[pop8Trans$transitions$fromClass == 2 & pop8Trans$transitions$toClass == 5 & pop8Trans$transitions$toColony == 1 ] 
  inf <- inf[!is.na(inf)]

  expect_equal( inf, 1 * 1 * 1 )

  # 2 -> 1 and 2 (i.e. class 3 to class 5)
  inf <- pop8Trans$transitions$rate[pop8Trans$transitions$fromClass == 3 & pop8Trans$transitions$toClass == 5 & pop8Trans$transitions$toColony == 1 ] 
  inf <- inf[!is.na(inf)]

  expect_equal( inf, 1 * 1 * 1 )

  
  # check coinfection  adds properly
  pop8 <- pop5
  pop8$I[c(2, 3), 1, 1] <- 1

  pop8Trans <- transRates(pop8, 1)

  inf <- pop8Trans$transitions$rate[pop8Trans$transitions$fromClass == 1 & pop8Trans$transitions$toClass == 2 & pop8Trans$transitions$toColony == 1 ] 
  inf <- inf[!is.na(inf)]


  expect_equal( inf, 1 * 1 * 1 )


})



test_that('initTransitions works.', {

  p <- makePop(model = 'SI', nPathogens = 2, nColonies = 2)

  p <- initTransitions(p)

  expect_true(all(is.na(p$transitions$rate[!p$transitions$type %in% c('death', 'birth')])))

  # transitions should be 1*2 birth, 2^n * 2 = 8 death, 4 * 2 infection, 2^n * 2 dispersal. 
  # Total = 26
  expect_equal(NROW(p$transitions), 26)
  
  expect_equal(p$transitions$type %>% table %>% as.vector, c(2, 8, 8, 8))

  # Just list the known answers for to and from

  trans1 <- p$transitions[p$transitions$fromColony == 1 & !is.na(p$transitions$fromColony), ] 


  # deaths from class 1-4, infection 1->2, 1->3, 2->4, 3->4, dispersal from class 1-4
  expect_equal(trans1$type, c(rep('death', 4), rep('infection', 4), rep('dispersal', 4)))

  expect_equal(trans1$fromClass, c(1:4, 1, 1, 2, 3, 1:4))

  # toColony. NA for death, 1 for infection, 2 for disperal
  expect_equal(trans1$toColony, c(rep(NA, 4), rep(1, 4), rep(2, 4)))

  # toClass. NA for death, 2, 3, 4, 4 (see above) then 1:4
  expect_equal(trans1$toClass, c(rep(NA, 4), 2,3,4,4, 1:4))



  # Test births seperately.
  expect_true(p$transitions[p$transitions$type == 'birth', 2:3] %>% is.na %>% all) 
  expect_equal(p$transitions[p$transitions$type == 'birth', 4:5] %>% unlist(use.names = FALSE), c(1, 2, 1, 1)) 
  

})



test_that('All rate calculating functions work.', {

  p <- makePop(model = 'SI', nPathogens = 2, nColonies = 2, birth = 1, transmission = 1, death = 1, dispersal = 1, crossImmunity = 1)

  # birthR

  expect_equal(birthR(p, 1), rep(10000, 2))

  # birth w/ 0 pop
  p$I[1, , 1] <- 0
  expect_equal(birthR(p, 1), rep(0, 2))

  # check adding
  p$I[, , 1] <- 1:8
  expect_equal(birthR(p, 1), c(sum(1:4), sum(5:8)))


  #deathR
  p$I[, , 1] <- 0
  expect_equal(deathR(p, 1), rep(0, 8))

  p$I[, , 1] <- 1:8 
  # fromColony is 1, 2, 1, 2. fromClass is 1,1,2,2
  #   so expect 1:8 by rows.
  expect_equal(deathR(p, 1), c(1,5,2,6,3,7,4,8) )

  # infectionR
  p$I[, , 1] <- 0
  expect_equal(infectionR(p, 1), rep(0, 4))

  p$I[1, , 1] <- 1:2
  expect_equal(infectionR(p, 1), rep(0, 4))

  p$I[1:2, , 1] <- c(0, 1, 0, 2) # zero susceptibles, infection in both colonies
  expect_equal(infectionR(p, 1), rep(0, 4))

  # both S and I classes positive.
  #   fromColony 1, 2, 1, 2. toClass 2, 2, 3, 3
  p$I[1:2, , 1] <- 1:4

  # So expect +ve infections first, then zero (second pathogen) infections after
  expect_equal(infectionR(p, 1), c(1*2, 3*4, 0, 0))

  
  # coinfectionR
  
  p$I[, , 1] <- 0
  expect_equal(coinfectionR(p, 1), rep(0, 4))

  p$I[2, , 1] <- 1:2
  expect_equal(coinfectionR(p, 1), rep(0, 4))

  p$I[2:3, , 1] <- 1:4
  expect_equal(coinfectionR(p, 1), c(1*2, 3*4, 1*2, 3*4))

  # check adding
  p$I[4, 1, 1] <- 1
  expect_equal(coinfectionR(p, 1), c(1 * (2 + 1), 3 * 4, (1 + 1) * 2, 3 * 4))

}) 



test_that('findInfectionTrans works', {

  expect_true(all.equal(infectionTrans(makePop(model = 'SI', nPathogens = 2)), cbind(c(1,1,2,3), c(2,3,4,4)), check.attributes=FALSE))

  expect_true(all.equal(infectionTrans(makePop(model = 'SI', nPathogens = 3)), cbind(c(1,1,1,2,3,2,4,3,4,5,6,7), c(2,3,4,5,5,6,6,7,7,8,8,8)), , check.attributes=FALSE))

})





test_that('findDisease added works', {

  # list of which disease is added in rows in p$transitions for coinfection
  expect_equal(findDiseaseAdded(makePop(model = 'SI', nPathogens = 2, nColonies = 2)), c(2,2,1,1))
  #expect_equal(findDiseaseAdded(makePop(model = 'SI', nPathogens = 3, nColonies = 2)), c(2,2,3,3,1,1,3,3,1,1,2,2,3,3,2,2,1,1))

})

test_that('Multiple infection death rate works', {

  # With issue #7, this used to crash.
  p <- makePop(events = 1e2+2, nColonies = 2, nPathogens = 2, model = 'SIS', infectDeath = 2, sample = 10)
  p <- seedPathogen(p, c(1:2), n = 700)
  p <- runSim(p)

  expect_true(exists('p'))
  


  
  # caused an error with old bug. #7 
  p <- makePop(events = 1e2, nColonies = 30, nPathogens = 3, model = 'SIS', recovery = 5,  sample = 10, dispersal = 0.05, birth = 0.05, death = 0.05, crossImmunity = 0.1, meanColonySize = 1000, infectDeath = 5)
  p <- seedPathogen(p, c(1:3), n = 700)

  expect_true(runSim(p), is_a("MetapopEpi"))


})











