context('Test makePop function that initializes populations.')

test_that('makePop creates correct data types', {
  p <- makePop()

  expect_equal(length(p), 16)
  expect_true(inherits(p$parameters, 'numeric'))
  expect_true(inherits(p$models, 'data.frame'))
  expect_true(inherits(p$I, 'array'))
  expect_true(inherits(p$edgeList, 'matrix'))
  expect_true(inherits(p$diseaseList, 'list'))
  expect_true(inherits(p$nClasses, 'integer'))
  expect_true(inherits(p$totalRate, 'numeric'))
  expect_true(length(p$totalRate) == 1)
  expect_true(inherits(p$randU, 'numeric'))
  expect_true(inherits(p$adjacency, 'matrix'))
  expect_true(inherits(p$diseaseClasses, 'character'))
  expect_true(inherits(p$whichClasses, 'matrix'))
  expect_true(inherits(p$transitions, 'data.frame'))
  expect_true(inherits(p$waiting, 'numeric'))
})

test_that('makePop creates correctly sized arrays for S and I', {
  p1 <- makePop(nColonies = 10, nPathogens = 5, events = 20)

  expect_true(all.equal(dim(p1$I), c(2^5, 10, 20 + 1)))

  # Currently nColonies = 1 crashes. Wasn't deliberate but is probably a good choice.
  expect_error(makePop(nColonies = 1, nPathogens = 5, events = 20))


  # nPathogens must be > 1 as well
  expect_error(makePop(nColonies = 10, nPathogens = 1, events = 20))


})

test_that('incorrect args gives errors', {
  
  expect_error(makePop(model = 'ERRRORR'))
  expect_error(makePop(model = 2))

  expect_error(makePop(nColonies = 1))
  expect_error(makePop(nColonies = 0))
  expect_error(makePop(nColonies = 4.3))

  expect_error(makePop(colonyDistr = 'ERRRORR'))

  expect_error(makePop(space = 0))
  expect_error(makePop(space = -1))
  expect_error(makePop(space = 'ERRRORR'))
  
  expect_error(makePop(maxDistance = 'ERRRORR'))
  expect_error(supressWarnings(makePop(maxDistance = 0)))
  expect_error(supressWarnings(makePop(maxDistance = -1)))

  expect_error(makePop(kernel = 'ERRRORR'))

  expect_error(makePop(events = 0))
  expect_error(makePop(events = -1))
  
  expect_error(makePop(colonySpatialDistr = 'ERRRORR'))

  expect_error(makePop(nPathogens = 0))

  expect_error(makePop(meanColonySize = -1))

  # 
  expect_error(makePop(birth = -1))
  expect_error(makePop(death = -1))
  expect_error(makePop(dispersal = -1))
  expect_error(makePop(transmission = -1))
  expect_error(makePop(recovery = -1))

  expect_error(makePop(crossImmunity = -0.1))
  expect_error(makePop(crossImmunity = 1.1))

})


