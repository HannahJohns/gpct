context("Calculation Tests for GPCT")

test_that("GPCT gives same results as genodds",{

  df <- alteplase
  df$treat <- as.numeric(df$treat)-1

  pwc <- gpct:::pairwise(sign(outer(df$mRS,df$mRS,"-")),attributes = data.frame(group = df$treat))

  z1 <- gpct:::gpct(pwc,explanatory_var = "group",splitTies = F)
  z2 <- genodds::genodds(df$mRS,df$treat,ties = "drop")

  expect_true(prod(abs(unlist(z1[1:3]) - c(z2$pooled_lnodds,z2$pooled_SElnodds, z2$pooled_p)) < 1e-6)==1)

})


test_that("GPCT gives the same results for summarised and non-summarised data",{

  df <- gpct::alteplase
  df$treat <- as.numeric(df$treat)-1

  pwc <- gpct:::pairwise(sign(outer(df$mRS,df$mRS,"-")),attributes = data.frame(group = df$treat))
  z1 <- gpct:::gpct(pwc,explanatory_var = "group",splitTies = F)

  df2 <- as.data.frame(table(df[,c("mRS","treat")]))
  df2$mRS <- as.numeric(as.character(df2$mRS))

  pwc2 <- gpct:::pairwise(sign(outer(df2$mRS,df2$mRS,"-")),
                          attributes = data.frame(group = df2$treat),
                          weight = df2$Freq)

  z2 <- gpct:::gpct(pwc2,explanatory_var="group",splitTies = F)

  expect_true(prod(abs(c(z1$logOR,z1$logSE,z1$pVal)-c(z2$logOR,z2$logSE,z2$pVal))<1e-6)==1)

})


test_that("GPCT returns the same results consistently",{

  censor <- matrix(0,3,3)
  censor[1,2] <- censor[1,3] <- 1


  # scales=NULL; rates=NULL; shapes=NULL; censor=NULL; heirarchy=NULL; treatRatio = 1
  # arrivalRate = NULL; maxt = NULL; missingAcrualRelax=NULL; minObsTime=0;
  # numThreads=NULL;
  # alpha=0.05; splitTies=FALSE;
  #

  scales =  rep(15,3)
  shapes =  rep(1,3)
  hr = rep(0.8,3)
  N=1200
  arrivalRate = 1.1
  maxt = 190
  missingAcrualRelax = "maxt"
  minObsTime = 10




  set.seed(2)
  gpct:::simulate.survival(scales = scales,
                         shapes = shapes,
                         censor = censor,
                         hr = hr,
                         N=N,
                         arrivalRate = arrivalRate,
                         maxt = maxt,
                         missingAcrualRelax = missingAcrualRelax,
                         minObsTime = minObsTime
  ) -> simulated_data

  repeat_calculation <- lapply(1:10,function(i){

    x <- gpct:::gpct(x = simulated_data$pwc,explanatory_var = "trt")

    c(x$logOR,x$logSE)

  })

  expect_true(do.call("all.equal",repeat_calculation))

})

test_that("GPCT and Win Ratio have practically the same effect and SE",{


  censor <- matrix(0,3,3)
  censor[1,2] <- censor[1,3] <- 1

  scales =  rep(15,3)
  shapes =  rep(1,3)
  hr = rep(0.8,3)
  N=1200
  arrivalRate = 1.1
  maxt = 190
  missingAcrualRelax = "maxt"
  minObsTime = 10


  set.seed(3)
  gpct:::simulate.survival(scales = scales,
                           shapes = shapes,
                           censor = censor,
                           hr = hr,
                           N=N,
                           arrivalRate = arrivalRate,
                           maxt = maxt,
                           missingAcrualRelax = missingAcrualRelax,
                           minObsTime = minObsTime
  ) -> simulated_data


  z1 <- gpct:::gpct(simulated_data$pwc,explanatory_var = "trt")
  z2 <- gpct:::winRatio(simulated_data$pwc,group = "trt")

  expect_true(prod(abs(unlist(z1[c("logOR","logSE")])-z2[[1]][c("logWR","seLogWR")])<1e-3)==1)

})

