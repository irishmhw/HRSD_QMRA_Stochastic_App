#' Annual Risk of Illness
#'
#' Calculate the annual risk of illness.
#' 
#' @param weightProbInf Daily Weighted Probability of Infection
#' @param illGivenInfection Pathogen specific risk of illness given infection.
#'
#' @return A vector of the annual risk of illness

AnnualRiskOfIllness <- function(weightProbInf, illGivenInfection){
  annRisk <- ifelse(weightProbInf < 10 ^ -6, 
                    weightProbInf * illGivenInfection * 365,
                    1 - ( 1 - weightProbInf * illGivenInfection) ^ 365)
  return(annRisk)
}

calcDALYs <- function(annRisk, DALYWeight){
  dalys <- annRisk * DALYWeight
  return(dalys)
}