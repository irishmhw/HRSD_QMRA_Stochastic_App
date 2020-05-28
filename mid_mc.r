source('./treatment_dist.r')
source('./gof.r')

qmra.exponential <- function(k,dose) {
  return(1 - exp(-k*dose))
}
qmra.bp <- function(alpha,N50,dose) {
  res <- 1 - ( 1 + ( dose / N50 ) * (2 ^ ( 1 / alpha) -1 ) ) ^ (-alpha)
  return(res)
}


qmra.MonteCarlo <- function(maxiter=100,
                            organisms=c('crypto','giardia','rota','campy','eColi'), # vector containing names
                            concentration,
                            coagName,
                            filtName,
                            disinfectName1,
                            disinfectName2
){
  
  treatment <- treatmentDist(maxiter)
  coag <- treatment$coag
  filtMap <- list(
    'filtNc' = treatment$filt_nc, 
    'sand' = treatment$filt_sand,
    'ultra' = treatment$filt_ultra,
    'micro' = treatment$filt_micro
  )
  filt <- filtMap[[filtName]]
  dInfMap <- list(
    'chlorine' = treatment$chlorine,
    'chloromines' = treatment$chloromines,
    'ozone' = treatment$ozone,
    'chlorine_dioxide' = treatment$chlorine_dioxide
  ) 
  disinfect1 <- dInfMap[[disinfectName1]]
  disinfect2 <- dInfMap[[disinfectName2]]
  #--- Build the matrices for the Monte Carlo ----
  rawConc <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  treatedConc <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  lr <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  dose <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  morbRatio <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  risk <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  riskMorb <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  
  ingestion <- matrix(nrow=maxiter, ncol=1)
  
  Crypto_k <- matrix(nrow=maxiter, ncol=1)
  Giardia_k <- matrix(nrow=maxiter, ncol=1)
  Rota_alpha <- matrix(nrow=maxiter, ncol=1)
  Rota_N50 <- matrix(nrow=maxiter, ncol=1)
  Campy_alpha <- matrix(nrow=maxiter, ncol=1)
  Campy_N50 <- matrix(nrow=maxiter, ncol=1)
  Ecoli_k <- matrix(nrow=maxiter, ncol=1)
  
  for(i in 1:maxiter){
    #Hard coding a bunch of values
    ingestion[i] <- TriRand(0.5, 1, 2.5)
    Crypto_k[i] <- TriRand(0.0074,0.0539,0.3044)                    # Allow uncertainty in k parameter
    morbRatio[,'crypto'][i] <- TriRand(0.3,0.5,0.7)                     # Uncertain morbidity ratio
    
    Giardia_k[i] <- TriRand(0.0085,0.0199,0.0371)                   # Allow uncertainty in k parameter
    morbRatio[,'giardia'][i] <- runif(1,0.4,0.4)                        # Uncertain morbidity ratio
    
    Rota_alpha[i] <- TriRand(1.64E-01,2.53E-01,5.18E-01)       # Allow uncertainty in alpha parameter
    Rota_N50[i] <- TriRand(2.49E+00,6.17E+00,1.89E+01)         # Allow uncertainty in alpha parameter		
    morbRatio[,'rota'][i] <- runif(1,0.88,0.88)                    # Uncertain morbidity ratio
    

    Campy_alpha[i] <- TriRand(4.99E-02,1.44E-01,2.66E-01)       # Allow uncertainty in alpha parameter
    Campy_N50[i] <- TriRand(8.11,8.9E+02,6.69E+03)         # Allow uncertainty in alpha parameter		
    morbRatio[,'campy'][i] <- runif(1,1,1)                    # Uncertain morbidity ratio
    
    Ecoli_k[i] <- TriRand(1.20E-04,2.18E-04,5.99E-04)       # Allow uncertainty in alpha parameter
    morbRatio[,'eColi'][i] <- runif(1,1,1)                    # Uncertain morbidity ratio
    
    for(org in organisms){
      rawConc[,org][i] <- concentration[,org][i]
      lr[,org][i] <- coag[,org][i] + filt[,org][i] + disinfect1[,org][i] + disinfect2[,org][i]
      treatedConc[,org][i] <- rawConc[,org][i] * 10 ^ (-1 * lr[,org][i])
      dose[,org][i] <- treatedConc[,org][i] * ingestion[i]
    }
    
    risk[,'crypto'][i] <- qmra.exponential(Crypto_k[i], dose[,'crypto'][i]) 
    risk[,'giardia'][i] <- qmra.exponential(Giardia_k[i], dose[,'giardia'][i])
    risk[,'eColi'][i] <- qmra.exponential(Ecoli_k[i], dose[,'eColi'][i])
    
    risk[,'rota'][i] <- qmra.bp(Rota_alpha[i], Rota_N50[i], dose[,'rota'][i])
    risk[,'campy'][i] <- qmra.bp(Campy_alpha[i], Campy_N50[i], dose[,'campy'][i])
    
    
    for(org in organisms){
      riskMorb[,org][i] <- risk[,org][i] * morbRatio[,org][i]
    }
  }
  
  return(list(dose=dose, risk=risk, lr=lr, treatedConc=treatedConc, riskMorb=riskMorb))
}
