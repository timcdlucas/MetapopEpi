context('All other exported functions.')


test_that('seedPathogen works correctly.', {
  set.seed(1)
  p <- makePop()

  s1 <- seedPathogen(p, 1)

  # Still should be 50000 individuals
  expect_true(sum(s1$I[,,1]) == 50000)
  
  # One colony should have lost a susceptible individual
  expect_true( sum(s1$I[1,,1] == 9999) == 1)
  
  # Same colony should have gained an infected individual
  lostSusc <- which(s1$I[1,,1] == 9999)
  gainedInf <- (which(s1$I[,,1] == 1, arr.ind = TRUE))[2]

  expect_equal(lostSusc, gainedInf)  

  # No coinfections should have occurred
  expect_true(all(s1$I[4:8,,] != 1))



  # Two pathogens
  # Doing pathogen 1 and path 3, just to check for weird stuff.
  s2 <- seedPathogen(p, c(1,3))

  # Still should be 50000 individuals
  expect_true(sum(s2$I[,,1]) == 50000)
  
  
  # Two susceptible individuals should have been lost
  # For some seeds, we get 1 x 9998, others 2 x 9999.
  # Would like to check former, but not sure how to force it. 
  expect_true( sum(s2$I[1,,1]) == 49998)

  expect_true( all(colSums(s2$I[,,1]) == 10000))
  
  # Same colony should have gained an infected individual
  lostSusc <- which(s2$I[1,,1] != 10000)
  gainedInf <- (which(s2$I[,,1] == 1, arr.ind = TRUE))[,2]

  # Possibly that both pathogens are seeded in same colony.
  # If so have to test differently.
  if(length(lostSusc) == length(gainedInf)){
    expect_equal(lostSusc, gainedInf)  
  } else {
    expect_true(all(lostSusc == gainedInf))
  }


  # No coinfections should have occurred
  expect_true(all(s2$I[5:8,,] != 1))

  # Check exactly one of each pathogen
  expect_true(sum(s2$I[2,,1]) == 1)
  expect_true(sum(s2$I[4,,1]) == 1)

  # All pathogens
  sAll <- seedPathogen(p, 1:p$parameters['nPathogens'])

  # Still should be 50000 individuals
  expect_true(sum(sAll$I[,,1]) == 50000)

  
  
  # Three susceptible individuals should have been lost
  # For some seeds, we get 1 x 9998, others 2 x 9999.
  # Would like to check former, but not sure how to force it. 
  expect_true( sum(sAll$I[1,,1]) == 49997)
  expect_true( all(colSums(sAll$I[,,1]) == 10000))
  

  # No coinfections should have occurred
  expect_true(all(s2$I[5:8,,] != 1))

  
  # Too many or too few pathogens to be seeded should cause error.
  expect_error(seedPathogen(p, 0))
  expect_error(seedPathogen(p, 4))
  

})


