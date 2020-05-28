# Header ------------------------------------------------------------------

#==============================================================================================================|
# Crypto_data_Dist_Fitting_GOF_MHW.R is an R (www.r-project.org) source code that operates a probability       |
# distribution optimization of data entered into R using the data vector. The optimzied probability            |
# distributions are then used to determine the best fitting models using the AIC weights, chosen to further    |
# discount the number of parameters in the probability distributions.                                          |
# The code is broken into sections due to length to allow for easier navigation                                |
#                                                                                                              |
# Coding and model redevelopment performed by Mark H. Weir Ph.D. of CAMRA Consultants LLC, NSF International   | 
# and College of Public Health and College of Engineering, The Ohio State University                           |
#                                                                                                              |
# All use and reproduction rights are reserved by Mark H. Weir Ph.D. and CAMRA Consultants LLC.                |
# CAMRA Consultants LLC. weirmarkh@gmail.com, camraconsultants@gmail.com; weir.95@osu.edu                      |
#==============================================================================================================|

# Load packages ---------------------------------------------------------
require(MASS)
# Load data and take log of data if needed and state distributions to be optimized ---------------
qmra.gof <- function(data, maxiter=10000){
  dists <- c("normal","lognormal","Weibull","geometric","exponential","logistic","Poisson","Cauchy")
  data <- data + 1e-11
  ldata <- log(data) + 1e-11

  fits = matrix(nrow=length(dists), ncol=1)

# Use fitdistr() to complete the probability distribution optimization ---------------------
for(i in 1:length(ldata)){
  SS_tot <- sum(data[i] - mean(data))
    
  normal <- fitdistr(data, "normal")
  RSS_norm <- (sum(data[i] - dnorm(data[i], normal$estimate[1], normal$estimate[2]))) ^ 2
	
  lognormal <- fitdistr(data, "lognormal")
  RSS_lnorm <- (sum(data[i] - dlnorm(data[i], lognormal$estimate[1], lognormal$estimate[2]))) ^ 2
	
  weibull <- fitdistr(data, "weibull")
  RSS_weibull <- (sum(data[i] - dweibull(data[i], weibull$estimate[1], weibull$estimate[2]))) ^ 2

  geometric <- fitdistr(data, "geometric")
  RSS_geometric <- (sum(data[i] - dgeom(data[i], geometric$estimate[1]))) ^ 2
	
  exponential <- fitdistr(data, "exponential")
  RSS_exponential <- (sum(data[i]- dexp(data[i], exponential$estimate[1]))) ^ 2
	
	logistic <- fitdistr(data, "logistic")
  RSS_logistic <- (sum(data[i] - dlogis(data[i], logistic$estimate[1], logistic$estimate[2]))) ^ 2
	
  poisson <- fitdistr(data, "poisson")
  RSS_poisson <- (sum(data[i] - dpois(data[i], poisson$estimate[1]))) ^ 2 

  cauchy <- fitdistr(data, "cauchy")
  RSS_cauchy <- (sum(data[i] - dcauchy(data[i], cauchy$estimate[1]))) ^ 2
	
}

# Perform a variety of goodness of fit and best fitting assessments AIC weights are used for decisions ----------

 RSS_all <- c( RSS_norm,  RSS_lnorm,  RSS_weibull,  RSS_geometric,  RSS_exponential,  RSS_logistic,  RSS_poisson,  RSS_cauchy)
 RSS_results <- data.frame(dists,  RSS_all)

 R_sqrd_norm <- SS_tot/ RSS_results[1,2];  
 R_sqrd_lnorm <- SS_tot/ RSS_results[2,2];  
 R_sqrd_Weibull <- SS_tot/ RSS_results[3,2];
 R_sqrd_geometric <- SS_tot/ RSS_results[4,2];  
 R_sqrd_exponential <- SS_tot/ RSS_results[5,2];  
 R_sqrd_logistic <- SS_tot/ RSS_results[6,2];
 R_sqrd_Poisson <- SS_tot/ RSS_results[7,2];  
 R_sqrd_Cauchy <- SS_tot/ RSS_results[8,2];
 Rsqrd_all <- c( R_sqrd_norm,  R_sqrd_lnorm,  R_sqrd_Weibull,  R_sqrd_geometric,  R_sqrd_exponential,  R_sqrd_logistic,  R_sqrd_Poisson,  R_sqrd_Cauchy)

 AIC_norm <- AIC( normal);  
 AIC_lnorm <- AIC( lognormal);  
 AIC_Weibull <- AIC( weibull);  
 AIC_geometric <- AIC( geometric);
 AIC_exponential <- AIC( exponential);  
 AIC_logistic <- AIC( logistic);  
 AIC_Poisson <- AIC( poisson);  
 AIC_Cauchy <- AIC( cauchy); 
 AIC_all <- c( AIC_norm,  AIC_lnorm,  AIC_Weibull,  AIC_geometric,  AIC_exponential,  AIC_logistic,  AIC_Poisson,  AIC_Cauchy)

 AICw <- exp((min( AIC_all)- AIC_all)/2)/sum(exp((min( AIC_all)- AIC_all)/2))

 BIC_norm <- BIC( normal);  
 BIC_lnorm <- BIC( lognormal);  
 BIC_Weibull <- BIC( weibull);  
 BIC_geometric <- BIC( geometric);
 BIC_exponential <- BIC( exponential);  
 BIC_logistic <- BIC( logistic);  
 BIC_Poisson <- BIC( poisson);  
 BIC_Cauchy <- BIC( cauchy); 
 BIC_all <- c( BIC_norm,  BIC_lnorm,  BIC_Weibull,  BIC_geometric,  BIC_exponential,  BIC_logistic,  BIC_Poisson,  BIC_Cauchy)

 Optim_Results <- data.frame(dists,  RSS_all,  Rsqrd_all,  AIC_all,  AICw,  BIC_all)
 
 
# colnames( Optim_Results) <- c("Distribution Name", "Residual Sum of Squares of Prediction", "R sqrd", "AIC", "AIC Weights", "BIC")
# write.csv( Optim_Results, file=" Sum_Sqrd_Error_Predctn_Raw_Data.csv") 


# Pull estimate values for inclusion in the Monte Carlo simulation ------------------------------------

estimate_1 <- c(normal$estimate[1], lognormal$estimate[1], weibull$estimate[1], geometric$estimate[1],
          exponential$estimate[1], logistic$estimate[1], poisson$estimate[1], cauchy$estimate[1])

 estimate_2 <- c(normal$estimate[2], lognormal$estimate[2], weibull$estimate[2], geometric$estimate[2],
                exponential$estimate[2], logistic$estimate[2], poisson$estimate[2], cauchy$estimate[2])

Results_for_MonteCarlo <- data.frame(dists,AIC_all, AICw, BIC_all, estimate_1, estimate_2)

Best_Fitting_Row <- which.max(Results_for_MonteCarlo$AICw)
Best_Fitting_Dist <- Results_for_MonteCarlo[Best_Fitting_Row,1]

# Model the concentrations within this code, until distname can be used in the Monte Carlo script ----------------
set.seed(37)
if(Best_Fitting_Dist=="normal"){
  distname <- rnorm; 
  Conc <- rnorm(maxiter,
    Results_for_MonteCarlo[Best_Fitting_Row,5], Results_for_MonteCarlo[Best_Fitting_Row,6]); 
  distdisplay <- "Normal"
}

if(Best_Fitting_Dist=="lognormal"){
  distname <- rlnorm;
  Conc <- rlnorm(maxiter,
    Results_for_MonteCarlo[Best_Fitting_Row,5], Results_for_MonteCarlo[Best_Fitting_Row,6]); distdisplay <- "Log Normal"
}

if(Best_Fitting_Dist=="Weibull"){
  distname <- rweibull;
  Conc <- rweibull(maxiter,
     Results_for_MonteCarlo[Best_Fitting_Row,5], Results_for_MonteCarlo[Best_Fitting_Row,6]); 
  distdisplay <- "Weibull"
}

if(Best_Fitting_Dist=="geometric"){
  distname <- rgeom; 
  Conc <- rgeom(maxiter,
     Results_for_MonteCarlo[Best_Fitting_Row,5]); 
  distdisplay <- "Geometric"
}

if(Best_Fitting_Dist=="exponential"){
  distname <- rexp;
  Conc <- rexp(maxiter,
     Results_for_MonteCarlo[Best_Fitting_Row,5]);  
  distdisplay <- "Exponential"
}

if(Best_Fitting_Dist=="logistic"){
  distname <- rlogis;
  Conc <- rlogis(maxiter,
     Results_for_MonteCarlo[Best_Fitting_Row,5], Results_for_MonteCarlo[Best_Fitting_Row,6]);  
  distdisplay <- "Logistic"
}

if(Best_Fitting_Dist=="Poisson"){
  distname <- rpois;
  Conc <- rpois(maxiter,
     Results_for_MonteCarlo[Best_Fitting_Row,5]);  
  distdisplay <- "Poisson"
}

if(Best_Fitting_Dist=="Cauchy"){
  distname <- rcauchy;
  Conc <- rcauchy(maxiter,
     Results_for_MonteCarlo[Best_Fitting_Row,5], Results_for_MonteCarlo[Best_Fitting_Row,6]);  
  distdisplay <- "Cauchy"
}
return(list(results=Results_for_MonteCarlo, conc=Conc, distDisplay=distdisplay ))
}