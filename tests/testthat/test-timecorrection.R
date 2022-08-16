context("Calculation Tests for simulating survival data")

test_that("Competing risks censoring behaves the same under transitive closure",{

  times <- matrix(rexp(100*3),100,3)

  censor <- matrix(0,3,3)
  censor[1,2] <- censor[2,3] <- 1

  correctedTimes1 <- matrix(0,100,3)
  gpct:::correctTimes(correctedTimes1,times,censor)

  censor[1,3] <- 1
  correctedTimes2 <- matrix(0,100,3)
  gpct:::correctTimes(correctedTimes2,times,censor)

  expect_equal(correctedTimes1,correctedTimes2)

})

