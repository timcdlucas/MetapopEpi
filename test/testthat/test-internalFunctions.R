context('Test internal functions for makePop and runnings sims.')


test_that('randEvent works correctly', {
 
  set.seed(2)
  pop <- makePop(nColonies = 2)
  pop <- seedPathogen(pop, c(1:3))
  
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
  pop <- makePop(events = 20)
  
  p1 <- runSim(pop)

  # All waiting times should be positive
  expect_true(all(p1$waiting[-1] > 0))


  # Check that the only thing changed is pop$waiting[t +1]
  pop2 <- makePop(events = 20)

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
  pop <- makePop(events = 20)

  pop2 <- transRates(pop, 1)

  expect_equal(pop$transitions, pop2$transitions)

  expect_true(all(pop2$transitions$rate[pop2$transitions$type == 'birth'] == 10))
  expect_true(all(pop2$transitions$rate[pop2$transitions$type == 'death' & pop2$transitions$fromClass == 1] == 10))
  expect_true(all(pop2$transitions$rate[pop2$transitions$type == 'death' & pop2$transitions$fromClass != 1] == 0))
  expect_true(all(pop2$transitions$rate[pop2$transitions$type == 'infection'] == 0))
  expect_true(all(pop2$transitions$rate[pop2$transitions$type == 'dispersal' & pop2$transitions$fromClass == 1] == 2.5))
  expect_true(all(pop2$transitions$rate[pop2$transitions$type == 'dispersal' & pop2$transitions$fromClass != 1] == 0))


  pop3 <- makePop(nColonies = 2, nPathogens = 3, transmission = 1, dispersal = 2, birth = 3, death = 4)

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
  

  pop5 <- makePop(nColonies = 2, nPathogens = 3, transmission = 1, crossImmunity = 1, meanColonySize = 1)

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
  pop8$I[c(1:3), 1, 1] <- c(0, 1, 1)x

  pop8Trans <- transRates(pop8, 1)

  pop8Trans$transitions$rate[pop8Trans$transitions$fromClass == 2 & pop8Trans$transitions$toClass == 6 & pop8Trans$transitions$toColony == 1 ] 
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




