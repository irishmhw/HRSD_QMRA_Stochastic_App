# This R source code performs a Monte Carlo written by Mark H. Weir PhD.
# This version of the code implements the TriRand function also written by 
# Mark H. Weir PhD.
require(ggplot2)
require(gridExtra)
require(reshape)
source("./TriRand.r")
source("./qmra_utils.r")
source('./treatment_dist.r')
source('./gof.r')

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
    'filt_nc' = treatment$filt_nc, 
    'sand' = treatment$filt_sand,
    'ultra' = treatment$filt_ultra,
    'micro' = treatment$filt_micro
  )
  filt <- filtMap[[filtName]]
  bioFilt <- treatment$bio_filt 
  
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
  
  ingestion <- rtri(maxiter, 0.5, 1, 2)
  
  Crypto_k <- rtri(maxiter, 0.0074, 0.0539, 0.3044)
  morbRatio[,'crypto'] <- rtri(maxiter, 0.3, 0.5, 0.7)
  
  Giardia_k <- rtri(maxiter, 0.0085,0.0199,0.0371)                   # Allow uncertainty in k parameter
  morbRatio[,'giardia'] <- runif(maxiter, 0.4, 0.4)                        # Uncertain morbidity ratio
  
  Rota_alpha <- rtri(maxiter, 1.64E-01, 2.53E-01, 5.18E-01)       # Allow uncertainty in alpha parameter
  Rota_N50 <- rtri(maxiter, 2.49E+00, 6.17E+00, 1.89E+01)         # Allow uncertainty in alpha parameter		
  morbRatio[,'rota']<- runif(maxiter,0.88,0.88)                    # Uncertain morbidity ratio
  
  Campy_alpha <- rtri(maxiter, 4.99E-02, 1.44E-01, 2.66E-01)       # Allow uncertainty in alpha parameter
  Campy_N50 <- rtri(maxiter, 8.11, 8.9E+02, 6.69E+03)         # Allow uncertainty in alpha parameter		
  morbRatio[,'campy'] <- runif(maxiter, 1, 1)                    # Uncertain morbidity ratio
  
  Ecoli_k <- rtri(maxiter, 1.20E-04, 2.18E-04, 5.99E-04)       # Allow uncertainty in alpha parameter
  morbRatio[,'eColi'] <- runif(maxiter, 1, 1)                    # Uncertain morbidity ratio
  
  
  for(org in organisms){
    coag[,org] <- coag[,org] * efficiency$coag
    disinfect1[,org] <- disinfect1[,org] * efficiency$disinfect1
    filt[,org] <- filt[,org] * efficiency$filt
    bioFilt[,org] <- bioFilt[,org] * efficiency$bioFilt
    disinfect2[,org] <- disinfect2[,org] * efficiency$disinfect2
    disinfect3[,org] <- disinfect3[,org] * efficiency$disinfect3
    
    rawConc[,org] <- concentration[,org]
    # present output for each lr
    # Order is coag, disinfect1, bioFilt, filt, disinfect2, disinfect3
    postCoagConc[,org] <- rawConc[,org] * 10 ^ (-1 * coag[,org])
    postDisinfect1Conc[,org] <- rawConc[,org] * 10 ^ (-1 * ( coag[,org] + disinfect1[,org]))
    postBioFiltConc[,org] <- rawConc[,org] * 10 ^ (-1 * ( coag[,org] + disinfect1[,org] + bioFilt[,org]))
    postFiltConc[,org] <- rawConc[,org] * 10 ^ (-1 * ( coag[,org] + disinfect1[,org] + bioFilt[,org] + filt[,org]))
    postDisinfect2Conc[,org] <- rawConc[,org] * 10 ^ (-1 * ( coag[,org] + disinfect1[,org] + bioFilt[,org] + filt[,org] + disinfect2[,org]))
    postDisinfect3Conc[,org] <- rawConc[,org] * 10 ^ (-1 * ( coag[,org] + disinfect1[,org] + bioFilt[,org] + filt[,org] + disinfect2[,org] +  disinfect3[,org]))
    
    postCoagDose[,org] <- postCoagConc[,org] * ingestion
    postDisinfect1Dose[,org] <- postDisinfect1Conc[,org] * ingestion
    postBioFiltDose[,org] <- postBioFiltConc[,org] * ingestion
    postFiltDose[,org] <- postFiltConc[,org] * ingestion
    postDisinfect2Dose[,org] <- postDisinfect2Conc[,org] * ingestion
    postDisinfect3Dose[,org] <- postDisinfect3Conc[,org] * ingestion
    
    treatedConc[,org] <- postDisinfect3Conc[,org] * ingestion
    dose[,org] <- treatedConc[,org] * ingestion
  }
  
  eps <- 1E-16
  # Coag
  postCoagRisk[,'crypto'] <- qmra.exponential(Crypto_k, postCoagDose[,'crypto']) + eps
  postCoagRisk[,'giardia'] <- qmra.exponential(Giardia_k, postCoagDose[,'giardia']) + eps
  postCoagRisk[,'eColi'] <- qmra.exponential(Ecoli_k, postCoagDose[,'eColi']) + eps
  
  postCoagRisk[,'rota'] <- qmra.bp(Rota_alpha, Rota_N50, postCoagDose[,'rota']) + eps
  postCoagRisk[,'campy'] <- qmra.bp(Campy_alpha, Campy_N50, postCoagDose[,'campy']) + eps
  
  # Disinfect1
  postDisinfect1Risk[,'crypto'] <- qmra.exponential(Crypto_k, postDisinfect1Dose[,'crypto']) + eps
  postDisinfect1Risk[,'giardia'] <- qmra.exponential(Giardia_k, postDisinfect1Dose[,'giardia']) + eps
  postDisinfect1Risk[,'eColi'] <- qmra.exponential(Ecoli_k, postDisinfect1Dose[,'eColi']) + eps
  
  postDisinfect1Risk[,'rota'] <- qmra.bp(Rota_alpha, Rota_N50, postDisinfect1Dose[,'rota']) + eps
  postDisinfect1Risk[,'campy'] <- qmra.bp(Campy_alpha, Campy_N50, postDisinfect1Dose[,'campy']) + eps
  
  # BioFilt
  postBioFiltRisk[,'crypto'] <- qmra.exponential(Crypto_k, postBioFiltDose[,'crypto']) + eps 
  postBioFiltRisk[,'giardia'] <- qmra.exponential(Giardia_k, postBioFiltDose[,'giardia']) + eps
  postBioFiltRisk[,'eColi'] <- qmra.exponential(Ecoli_k, postBioFiltDose[,'eColi']) + eps
  
  postBioFiltRisk[,'rota'] <- qmra.bp(Rota_alpha, Rota_N50, postBioFiltDose[,'rota']) + eps
  postBioFiltRisk[,'campy'] <- qmra.bp(Campy_alpha, Campy_N50, postBioFiltDose[,'campy']) + eps
  
  # Filt
  postFiltRisk[,'crypto'] <- qmra.exponential(Crypto_k, postFiltDose[,'crypto']) + eps 
  postFiltRisk[,'giardia'] <- qmra.exponential(Giardia_k, postFiltDose[,'giardia']) + eps
  postFiltRisk[,'eColi'] <- qmra.exponential(Ecoli_k, postFiltDose[,'eColi']) + eps
  
  postFiltRisk[,'rota'] <- qmra.bp(Rota_alpha, Rota_N50, postFiltDose[,'rota']) + eps
  postFiltRisk[,'campy'] <- qmra.bp(Campy_alpha, Campy_N50, postFiltDose[,'campy']) + eps
  
  
  # Disinfect2
  postDisinfect2Risk[,'crypto'] <- qmra.exponential(Crypto_k, postDisinfect2Dose[,'crypto']) + eps 
  postDisinfect2Risk[,'giardia'] <- qmra.exponential(Giardia_k, postDisinfect2Dose[,'giardia']) + eps
  postDisinfect2Risk[,'eColi'] <- qmra.exponential(Ecoli_k, postDisinfect2Dose[,'eColi']) + eps
  
  postDisinfect2Risk[,'rota'] <- qmra.bp(Rota_alpha, Rota_N50, postDisinfect2Dose[,'rota']) + eps
  postDisinfect2Risk[,'campy'] <- qmra.bp(Campy_alpha, Campy_N50, postDisinfect2Dose[,'campy']) + eps
  
  # Disinfect3
  postDisinfect3Risk[,'crypto'] <- qmra.exponential(Crypto_k, postDisinfect3Dose[,'crypto']) + eps 
  postDisinfect3Risk[,'giardia'] <- qmra.exponential(Giardia_k, postDisinfect3Dose[,'giardia']) + eps
  postDisinfect3Risk[,'eColi'] <- qmra.exponential(Ecoli_k, postDisinfect3Dose[,'eColi']) + eps
  
  postDisinfect3Risk[,'rota'] <- qmra.bp(Rota_alpha, Rota_N50, postDisinfect3Dose[,'rota']) + eps
  postDisinfect3Risk[,'campy'] <- qmra.bp(Campy_alpha, Campy_N50, postDisinfect3Dose[,'campy']) + eps
  
  # final
  risk[,'crypto'] <- qmra.exponential(Crypto_k, dose[,'crypto']) + eps 
  risk[,'giardia'] <- qmra.exponential(Giardia_k, dose[,'giardia']) + eps
  risk[,'eColi'] <- qmra.exponential(Ecoli_k, dose[,'eColi']) + eps
  
  risk[,'rota'] <- qmra.bp(Rota_alpha, Rota_N50, dose[,'rota']) + eps
  risk[,'campy'] <- qmra.bp(Campy_alpha, Campy_N50, dose[,'campy']) + eps
  
  
  for(org in organisms){
    postCoagMorb[,org] <- postCoagRisk[,org] * morbRatio[,org]
    postDisinfect1Morb[,org] <- postDisinfect1Risk[,org] * morbRatio[,org]
    postBioFiltMorb[,org] <- postBioFiltRisk[,org] * morbRatio[,org]
    postFiltMorb[,org] <- postFiltRisk[,org] * morbRatio[,org]
    postDisinfect2Morb[,org] <- postDisinfect2Risk[,org] * morbRatio[,org]
    postDisinfect3Morb[,org] <- postDisinfect2Risk[,org] * morbRatio[,org]
    
    riskMorb[,org] <- risk[,org] * morbRatio[,org]
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
      fun.y = median,
      geom = "point",
      size = 2.5,
      color = "white"
    ) +
    stat_summary(
      fun.y = mean,
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
      fun.y = median,
      geom = "point",
      size = 2.5,
      color = "white"
    ) +
    stat_summary(
      fun.y = mean,
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
      fun.y = median,
      geom = "point",
      size = 2.5,
      color = "white"
    ) +
    stat_summary(
      fun.y = mean,
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
      fun.y = median,
      geom = "point",
      size = 2.5,
      color = "white"
    ) +
    stat_summary(
      fun.y = mean,
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
      fun.y = median,
      geom = "point",
      size = 2.5,
      color = "white"
    ) +
    stat_summary(
      fun.y = mean,
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
      fun.y = median,
      geom = "point",
      size = 2.5,
      color = "white"
    ) +
    stat_summary(
      fun.y = mean,
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