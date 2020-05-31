library(testthat)

source('./monte_carlo.r')

maxiter <- 1000 # should be 1000

# use default data for most tests
defaultData <- read.csv('default_MPN.csv', header=TRUE, col.names = c('x'))

cryptoData <- defaultData
giardiaData <- defaultData
rotaData <- defaultData
campyData <- defaultData
eColiData <- defaultData

getDataInput <- function(){
  return(data.frame(
    crypto = cryptoData$x,
    giardia = giardiaData$x,
    rota = rotaData$x,
    campy = campyData$x,
    eColi = eColiData$x
  ))
}

d <- getDataInput()

res <- qmra.MonteCarlo(d, 
                       maxiter=maxiter, 
                       filtName='filt_nc',
                       disinfectName1='ozone', 
                       disinfectName2 = 'uv', 
                       disinfectName3 = 'chlorine', 
                       ingest = list(min=0.5, mean=1, max=2.5),
                       efficiency = list(
                         coag=1,
                         disinfect1=1,
                         filt=1,
                         bioFilt=1,
                         disinfect2=1,
                         disinfect3=1
                       )
                       )

# Test a number of pathogens at diffent stages.
prbs <- c(0.005, 0.5, 0.995)

test_that('Post coag illness quantiles are in expected range', {
  
  pcg <- apply(res$postCoagMorb, 2, quantile, probs=prbs)
  
  # Crypto PostCoag
  expect_equal(pcg[1,1], 2.12e-6, tolerance=1e-6)
  expect_equal(pcg[2,1], 9.19e-4, tolerance=1e-4)
  expect_equal(pcg[3,1], 8.55e-2, tolerance=1e-2)
  
  
  # Ecoli PostCoag
  ecoli <- c(3.47e-09, 2.98e-06, 2.60e-04)
  expect_equal(pcg[1,5], ecoli[1], tolerance=1e-8)
  expect_equal(pcg[2,5], ecoli[2], tolerance=1e-5)
  expect_equal(pcg[3,5], ecoli[3], tolerance=1e-3)
  
})

test_that('...post ozone illness...',{
  
  poz <- apply(res$postDisinfect1Risk, 2, quantile, probs=prbs)
  
  crypto <- c(2.293685e-07, 9.951544e-05, 9.970569e-03)
  expect_equal(poz[1,1], crypto[1], tolerance=1e-7)
  expect_equal(poz[2,1], crypto[2], tolerance=1e-5)
  expect_equal(poz[3,1], crypto[3], tolerance=1e-2)
})

test_that('...post biofilt illness...',{
  
  pbf <- apply(res$postBioFiltMorb, 2, quantile, probs=prbs)
  
  giardia <- c(2.211e-09, 2.89e-07, 3.15e-06)
  expect_equal(pbf[1,2], giardia[1], tolerance=1e-9)
  expect_equal(pbf[2,2], giardia[2], tolerance=1e-7)
  expect_equal(pbf[3,2], giardia[3], tolerance=1e-6)
  
  rota <- c(1.74e-10, 1.28e-07, 2.13e-05 )
  expect_equal(pbf[1,3], rota[1], tolerance=1e-9)
  expect_equal(pbf[2,3], rota[2], tolerance=1e-7)
  expect_equal(pbf[3,3], rota[3], tolerance=1e-5)
  
})

test_that('Post UV illness quantiles are in expected range' {
  
  puv <- apply(res$postDisinfect2Morb, 2, quantile, probs=prbs)

  # Rota PostUV
  expect_equal(puv[1,3], 2.34e-14, tolerance=1e-14)
  expect_equal(puv[2,3], 2.10e-11, tolerance=1e-10)
  expect_equal(puv[3,3], 4.92e-09, tolerance=1e-08)
  
  # Giardia PostUV
  giardia <- c(7.37e-14, 1.85e-11, 7.21e-10)
  expect_equal(puv[1,4], giardia[1], tolerance=1e-14)
  expect_equal(puv[2,4], giardia[2], tolerance=1e-10)
  expect_equal(puv[3,4], giardia[3], tolerance=1e-08)
  
})

test_that('Final illness quantiles are in expected range' {
  
  rqs <- apply(res$risk, 2, quantile, probs=prbs)
  
  # Rota PostUV
  expect_equal(rqs[1,1],2.69e-08)
  expect_equal(rqs[2,3],2.10e-11, tolerance=1e-10)
  expect_equal(rqs[3,3],4.92e-09, tolerance=1e-8)
})
