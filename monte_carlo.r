# This R source code performs a Monte Carlo written by Mark H. Weir PhD.
# This version of the code implements the TriRand function also written by 
# Mark H. Weir PhD.
require(ggplot2)
require(gridExtra)
require(reshape)
source("./TriRand.r")
source("./qmra_utils.r")

qmra.MonteCarlo <- function(inputData, 
                               maxiter = 1000,
                               coagName = 'coag',
                               filtName = 'sand',
                               disinfectName1 = 'chlorine',
                               disinfectName2 = 'ozone',
                               disinfectName3 = 'uv',
                               organisms = c('crypto','giardia','rota','campy','eColi'),
                               efficiency,
                               ingest
                            ){
  concentration <- data.frame(
    crypto=rep(NA, maxiter),
    giardia=rep(NA, maxiter),
    rota=rep(NA, maxiter),
    campy=rep(NA, maxiter),
    eColi=rep(NA, maxiter)
  )
  for (org in organisms){
    gofRes <- suppressWarnings(qmra.gof(inputData[,org], maxiter=maxiter))
    concentration[,org] <- gofRes$conc
  }
  res <- qmra._MonteCarlo(maxiter, 
                         concentration = concentration, 
                         coagName = coagName, 
                         filtName = filtName, 
                         disinfectName1 = disinfectName1, 
                         disinfectName2 = disinfectName2,
                         disinfectName3 = disinfectName3,
                         efficiency = efficiency,
                         ingest = ingest
          )
  
  
  return(res)
}

# _MonteCarlo
#
# Protected function where most of the monte carlo is performed
#

qmra._MonteCarlo <- function(maxiter=100,
                            organisms=c('crypto','giardia','rota','campy','eColi'), # vector containing names
                            concentration,
                            coagName,
                            filtName,
                            disinfectName1,
                            disinfectName2,
                            disinfectName3,
                            efficiency,
                            ingest
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
  bioFilt <- treatment$bio_filt #FIXME ******
  
  dInfMap <- list(
    'chlorine' = treatment$chlorine,
    'chloromines' = treatment$chloromines,
    'ozone' = treatment$ozone,
    'chlorine_dioxide' = treatment$chlorine_dioxide,
    'uv' = treatment$uv
  ) 
  
  disinfect1 <- dInfMap[[disinfectName1]]
  disinfect2 <- dInfMap[[disinfectName2]]
  disinfect3 <- dInfMap[[disinfectName3]]
  
  #--- Build the matrices for the Monte Carlo ----
  rawConc <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postCoagConc <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postFiltConc <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postBioFiltConc <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postDisinfect1Conc <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postDisinfect2Conc <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postDisinfect3Conc <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  treatedConc <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  lr <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  
  
  postCoagDose <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postFiltDose <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postBioFiltDose <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postDisinfect1Dose <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postDisinfect2Dose <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postDisinfect3Dose <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  dose <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  
  
  morbRatio <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  
  postCoagRisk <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postFiltRisk <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postBioFiltRisk <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postDisinfect1Risk <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postDisinfect2Risk <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postDisinfect3Risk <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  risk <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  
  postCoagMorb <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postFiltMorb <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postBioFiltMorb <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postDisinfect1Morb <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postDisinfect2Morb <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
  postDisinfect3Morb <- matrix(nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
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
    
    ingestion[i] <- TriRand(ingest$min, ingest$mean, ingest$max)
    
    #Hard coding a bunch of values
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
      coag[,org][i] <- coag[,org][i] * efficiency$coag
      disinfect1[,org][i] <- disinfect1[,org][i] * efficiency$disinfect1
      filt[,org][i] <- filt[,org][i] * efficiency$filt
      bioFilt[,org][i] <- bioFilt[,org][i] * efficiency$bioFilt
      disinfect2[,org][i] <- disinfect2[,org][i] * efficiency$disinfect2
      disinfect3[,org][i] <- disinfect3[,org][i] * efficiency$disinfect3
      
      rawConc[,org][i] <- concentration[,org][i]
      # present output for each lr
      # Order is coag, disinfect1, bioFilt, filt, disinfect2, disinfect3
      postCoagConc[,org][i] <- rawConc[,org][i] * 10 ^ (-1 * coag[,org][i])
      postDisinfect1Conc[,org][i] <- rawConc[,org][i] * 10 ^ (-1 * ( coag[,org][i] + disinfect1[,org][i]))
      postBioFiltConc[,org][i] <- rawConc[,org][i] * 10 ^ (-1 * ( coag[,org][i] + disinfect1[,org][i] + bioFilt[,org][i]))
      postFiltConc[,org][i] <- rawConc[,org][i] * 10 ^ (-1 * ( coag[,org][i] + disinfect1[,org][i] + bioFilt[,org][i] + filt[,org][i]))
      postDisinfect2Conc[,org][i] <- rawConc[,org][i] * 10 ^ (-1 * ( coag[,org][i] + disinfect1[,org][i] + bioFilt[,org][i] + filt[,org][i] + disinfect2[,org][i]))
      postDisinfect3Conc[,org][i] <- rawConc[,org][i] * 10 ^ (-1 * ( coag[,org][i] + disinfect1[,org][i] + bioFilt[,org][i] + filt[,org][i] + disinfect2[,org][i] +  disinfect3[,org][i]))
      
      postCoagDose[,org][i] <- postCoagConc[,org][i] * ingestion[i]
      postDisinfect1Dose[,org][i] <- postDisinfect1Conc[,org][i] * ingestion[i]
      postBioFiltDose[,org][i] <- postBioFiltConc[,org][i] * ingestion[i]
      postFiltDose[,org][i] <- postFiltConc[,org][i] * ingestion[i]
      postDisinfect2Dose[,org][i] <- postDisinfect2Conc[,org][i] * ingestion[i]
      postDisinfect3Dose[,org][i] <- postDisinfect3Conc[,org][i] * ingestion[i]
      
      treatedConc[,org][i] <- postDisinfect3Conc[,org][i] * ingestion[i]
      dose[,org][i] <- treatedConc[,org][i] * ingestion[i]
    }
    
    eps <- 1E-16
    # Coag
    postCoagRisk[,'crypto'][i] <- qmra.exponential(Crypto_k[i], postCoagDose[,'crypto'][i]) + eps
    postCoagRisk[,'giardia'][i] <- qmra.exponential(Giardia_k[i], postCoagDose[,'giardia'][i]) + eps
    postCoagRisk[,'eColi'][i] <- qmra.exponential(Ecoli_k[i], postCoagDose[,'eColi'][i]) + eps
    
    postCoagRisk[,'rota'][i] <- qmra.bp(Rota_alpha[i], Rota_N50[i], postCoagDose[,'rota'][i]) + eps
    postCoagRisk[,'campy'][i] <- qmra.bp(Campy_alpha[i], Campy_N50[i], postCoagDose[,'campy'][i]) + eps
    
    # Disinfect1
    postDisinfect1Risk[,'crypto'][i] <- qmra.exponential(Crypto_k[i], postDisinfect1Dose[,'crypto'][i]) + eps
    postDisinfect1Risk[,'giardia'][i] <- qmra.exponential(Giardia_k[i], postDisinfect1Dose[,'giardia'][i]) + eps
    postDisinfect1Risk[,'eColi'][i] <- qmra.exponential(Ecoli_k[i], postDisinfect1Dose[,'eColi'][i]) + eps
    
    postDisinfect1Risk[,'rota'][i] <- qmra.bp(Rota_alpha[i], Rota_N50[i], postDisinfect1Dose[,'rota'][i]) + eps
    postDisinfect1Risk[,'campy'][i] <- qmra.bp(Campy_alpha[i], Campy_N50[i], postDisinfect1Dose[,'campy'][i]) + eps
    
    # BioFilt
    postBioFiltRisk[,'crypto'][i] <- qmra.exponential(Crypto_k[i], postBioFiltDose[,'crypto'][i]) + eps 
    postBioFiltRisk[,'giardia'][i] <- qmra.exponential(Giardia_k[i], postBioFiltDose[,'giardia'][i]) + eps
    postBioFiltRisk[,'eColi'][i] <- qmra.exponential(Ecoli_k[i], postBioFiltDose[,'eColi'][i]) + eps
    
    postBioFiltRisk[,'rota'][i] <- qmra.bp(Rota_alpha[i], Rota_N50[i], postBioFiltDose[,'rota'][i]) + eps
    postBioFiltRisk[,'campy'][i] <- qmra.bp(Campy_alpha[i], Campy_N50[i], postBioFiltDose[,'campy'][i]) + eps
    
    # Filt
    postFiltRisk[,'crypto'][i] <- qmra.exponential(Crypto_k[i], postFiltDose[,'crypto'][i]) + eps 
    postFiltRisk[,'giardia'][i] <- qmra.exponential(Giardia_k[i], postFiltDose[,'giardia'][i]) + eps
    postFiltRisk[,'eColi'][i] <- qmra.exponential(Ecoli_k[i], postFiltDose[,'eColi'][i]) + eps
    
    postFiltRisk[,'rota'][i] <- qmra.bp(Rota_alpha[i], Rota_N50[i], postFiltDose[,'rota'][i]) + eps
    postFiltRisk[,'campy'][i] <- qmra.bp(Campy_alpha[i], Campy_N50[i], postFiltDose[,'campy'][i]) + eps
    
    
    # Disinfect2
    postDisinfect2Risk[,'crypto'][i] <- qmra.exponential(Crypto_k[i], postDisinfect2Dose[,'crypto'][i]) + eps 
    postDisinfect2Risk[,'giardia'][i] <- qmra.exponential(Giardia_k[i], postDisinfect2Dose[,'giardia'][i]) + eps
    postDisinfect2Risk[,'eColi'][i] <- qmra.exponential(Ecoli_k[i], postDisinfect2Dose[,'eColi'][i]) + eps
    
    postDisinfect2Risk[,'rota'][i] <- qmra.bp(Rota_alpha[i], Rota_N50[i], postDisinfect2Dose[,'rota'][i]) + eps
    postDisinfect2Risk[,'campy'][i] <- qmra.bp(Campy_alpha[i], Campy_N50[i], postDisinfect2Dose[,'campy'][i]) + eps
    
    # Disinfect3
    postDisinfect3Risk[,'crypto'][i] <- qmra.exponential(Crypto_k[i], postDisinfect3Dose[,'crypto'][i]) + eps 
    postDisinfect3Risk[,'giardia'][i] <- qmra.exponential(Giardia_k[i], postDisinfect3Dose[,'giardia'][i]) + eps
    postDisinfect3Risk[,'eColi'][i] <- qmra.exponential(Ecoli_k[i], postDisinfect3Dose[,'eColi'][i]) + eps
    
    postDisinfect3Risk[,'rota'][i] <- qmra.bp(Rota_alpha[i], Rota_N50[i], postDisinfect3Dose[,'rota'][i]) + eps
    postDisinfect3Risk[,'campy'][i] <- qmra.bp(Campy_alpha[i], Campy_N50[i], postDisinfect3Dose[,'campy'][i]) + eps
    
    # final
    risk[,'crypto'][i] <- qmra.exponential(Crypto_k[i], dose[,'crypto'][i]) + eps 
    risk[,'giardia'][i] <- qmra.exponential(Giardia_k[i], dose[,'giardia'][i]) + eps
    risk[,'eColi'][i] <- qmra.exponential(Ecoli_k[i], dose[,'eColi'][i]) + eps
    
    risk[,'rota'][i] <- qmra.bp(Rota_alpha[i], Rota_N50[i], dose[,'rota'][i]) + eps
    risk[,'campy'][i] <- qmra.bp(Campy_alpha[i], Campy_N50[i], dose[,'campy'][i]) + eps
    
    
    for(org in organisms){
      postCoagMorb[,org][i] <- postCoagRisk[,org][i] * morbRatio[,org][i]
      postDisinfect1Morb[,org][i] <- postDisinfect1Risk[,org][i] * morbRatio[,org][i]
      postBioFiltMorb[,org][i] <- postBioFiltRisk[,org][i] * morbRatio[,org][i]
      postFiltMorb[,org][i] <- postFiltRisk[,org][i] * morbRatio[,org][i]
      postDisinfect2Morb[,org][i] <- postDisinfect2Risk[,org][i] * morbRatio[,org][i]
      postDisinfect3Morb[,org][i] <- postDisinfect2Risk[,org][i] * morbRatio[,org][i]
      
      riskMorb[,org][i] <- risk[,org][i] * morbRatio[,org][i]
    }
  }
  
  return(list(
      postCoagDose=postCoagDose,
      postBioFiltDose=postBioFiltDose,
      postFiltDose=postFiltDose,
      postDisinfect1Dose=postDisinfect1Dose,
      postDisinfect2Dose=postDisinfect2Dose,
      postDisinfect3Dose=postDisinfect3Dose,
      dose=dose,
      postCoagRisk=postCoagRisk,
      postBioFiltRisk=postBioFiltRisk,
      postFiltRisk=postFiltRisk,
      postDisinfect1Risk=postDisinfect1Risk,
      postDisinfect2Risk=postDisinfect2Risk,
      postDisinfect3Risk=postDisinfect3Risk,
      risk=risk, 
      lr=lr, 
      treatedConc=treatedConc, 
      postCoagMorb=postCoagMorb,
      postBioFiltMorb=postBioFiltMorb,
      postFiltMorb=postFiltMorb,
      postDisinfect1Morb=postDisinfect1Morb,
      postDisinfect2Morb=postDisinfect2Morb,
      postDisinfect3Morb=postDisinfect3Morb,
      riskMorb=riskMorb
    ))
}


# MonteCarloPlot
# 
# Produces histograms and scatterplots based on result object from qmra.MonteCarlo.
#
qmra.MonteCarloPlot <- function(mcResults,
                    organisms = c('crypto','giardia','rota','campy','eColi'),
                    orgColor = list(
                      crypto = 'red',
                      giardia = 'sky blue',
                      rota = 'spring green',
                      campy = 'orange',
                      eColi = 'blue violet'
                    ),
                    nbins = 15
                    ){
  
  dose <- data.frame(mcResults$dose)
  risk <- data.frame(mcResults$risk)
  riskMorb <- data.frame(mcResults$riskMorb)
  
  plots <- list()
  i <- 1
  for (org in organisms){
    hist1 <- ggplot(risk,aes(x='risk')) +
      geom_histogram(aes(x=risk[,org]), fill=orgColor[[org]], bins=nbins) +
      geom_rug(aes(x=risk[,org]), color=orgColor[[org]]) + 
      labs(x='Risk')
    
    point1 <- ggplot(risk, aes(x=risk[,org])) + 
      stat_bin(aes(y = (..count..)/sum(..count..)), bins=nbins, geom='point',color=orgColor[[org]]) +
      labs(x='Risk', y='Percent')
    
    hist2 <- ggplot(riskMorb,aes(x='riskMorb')) +
      geom_histogram(aes(x=riskMorb[,org]),fill=orgColor[[org]], bins=nbins) +
      geom_rug(aes(x=riskMorb[,org]), color=orgColor[[org]]) + 
      labs(x='Risk Morbidity')
    
    point2 <- ggplot(risk, aes(x=riskMorb[,org])) + 
      stat_bin(aes(y = (..count..)/sum(..count..)), bins=nbins, geom='point',color=orgColor[[org]]) +
      labs(x='Risk Morbidity', y='Percent')
    
    plots[[i]] <- arrangeGrob(hist1, hist2, point1, point2, nrow=2, ncol=2, top=org)
    i <- i + 1
  }
  
  grid.arrange(grobs=plots, nrows=5, ncol=1, heights=unit(rep(4,5), rep('in', 5)))
}


qmra.treatmentPlots <- function(orgName, doses, risks, illness){
  
  # BOXPLOT and VIOLIN PLOT
  doses_ggplot <- suppressMessages(melt(doses))
  colnames(doses_ggplot) <- c("Treatment", "Dose")
  
  doses_box <- ggplot(doses_ggplot, aes(Treatment, Dose, fill = Treatment)) +
    geom_boxplot() +
    scale_y_log10() +
    theme(axis.text.x = element_text(
      face = "bold",
      color = "black",
      size = 10,
      angle = 45
    )) +
    labs(x = c(expression(bold(
      "Treatement Process" ^ a
    ))),
    y = expression(bold( ~ log[10] ~ of ~ Exiting ~ Dose))) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 2.5,
      color = "white"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 2,
      color = "black",
      shape = 15
    ) +
    guides(fill = guide_legend(override.aes = list(shape = NA)))
  
  doses_violin <-
    ggplot(doses_ggplot, aes(Treatment, Dose, fill = Treatment)) +
    geom_violin() +
    scale_y_log10() +
    theme(axis.text.x = element_text(
      face = "bold",
      color = "black",
      size = 10,
      angle = 45
    )) +
    labs(x = c(expression(bold(
      "Treatement Process" ^ a
    ))),
    y = expression(bold( ~ log[10] ~ of ~ Exiting ~ Dose))) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 2.5,
      color = "white"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 2,
      color = "black",
      shape = 15
    ) +
    guides(fill = guide_legend(override.aes = list(shape = NA)))
  
  risks_ggplot <-
    suppressMessages(melt(risks))
  colnames(risks_ggplot) <- c("Treatment", "Risk")
  
  risks_box <-
    ggplot(risks_ggplot, aes(Treatment, Risk, fill = Treatment)) +
    geom_boxplot() +
    scale_y_log10() +
    theme(axis.text.x = element_text(
      face = "bold",
      color = "black",
      size = 10,
      angle = 45
    )) +
    labs(x = c(expression(bold(
      "Treatement Process" ^ a
    ))),
    y = expression(bold( ~ log[10] ~ of ~ Exiting ~ Risk ~ of ~
                           Infection))) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 2.5,
      color = "white"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 2,
      color = "black",
      shape = 15
    ) +
    guides(fill = guide_legend(override.aes = list(shape = NA)))
  
  risks_violin <-
    ggplot(risks_ggplot, aes(Treatment, Risk, fill = Treatment)) +
    geom_violin() +
    scale_y_log10() +
    theme(axis.text.x = element_text(
      face = "bold",
      color = "black",
      size = 10,
      angle = 45
    )) +
    labs(x = c(expression(bold(
      "Treatement Process" ^ a
    ))),
    y = expression(bold( ~ log[10] ~ of ~ Exiting ~ Risk ~ of ~
                           Infection))) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 2.5,
      color = "white"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 2,
      color = "black",
      shape = 15
    ) +
    guides(fill = guide_legend(override.aes = list(shape = NA)))
  
  illness_ggplot <- suppressMessages(melt(illness))
  colnames(illness_ggplot) <- c("Treatment", "Illness")
  
  illness_box <- ggplot(illness_ggplot, aes(Treatment, Illness, fill = Treatment)) +
    geom_boxplot() +
    scale_y_log10() +
    theme(axis.text.x = element_text(
      face = "bold",
      color = "black",
      size = 10,
      angle = 45
    )) +
    labs(x = c(expression(bold(
      "Treatement Process" ^ a
    ))),
    y = expression(bold( ~ log[10] ~ of ~ Exiting ~ Risk ~ of ~
                           Illness))) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 2.5,
      color = "white"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 2,
      color = "black",
      shape = 15
    ) +
    guides(fill = guide_legend(override.aes = list(shape = NA)))
  
  illness_violin <-
    ggplot(illness_ggplot, aes(Treatment, Illness, fill = Treatment)) +
    geom_violin() +
    scale_y_log10() +
    theme(axis.text.x = element_text(
      face = "bold",
      color = "black",
      size = 10,
      angle = 45
    )) +
    labs(x = c(expression(bold(
      "Treatement Process" ^ a
    ))),
    y = expression(bold( ~ log[10] ~ of ~ Exiting ~ Risk ~ of ~
                           Illness))) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 2.5,
      color = "white"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 2,
      color = "black",
      shape = 15
    ) +
    guides(fill = guide_legend(override.aes = list(shape = NA)))
  
  plots <- arrangeGrob(
    #doses_violin, 
    #doses_box, 
    risks_violin,
    risks_box,
    illness_violin,
    illness_box, nrow=2, ncol=2, top=orgName)
  
  return(plots)
}