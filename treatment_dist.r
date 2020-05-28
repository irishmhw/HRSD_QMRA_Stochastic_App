
# Header ------------------------------------------------------------------

#==============================================================================================================|
# Treatment_Distributions.R is an R (www.r-project.org) source code that randomly samples from distributions   |
# for specified treatment processes and pathogens.                                                             |
# The code is broken into sections due to length to allow for easier navigation                                |
#                                                                                                              |
# Coding and model redevelopment performed by Mark H. Weir Ph.D. of CAMRA Consultants LLC, NSF International   | 
# and College of Public Health and College of Engineering, The Ohio State University                           |
#                                                                                                              |
# All use and reproduction rights are reserved by Mark H. Weir Ph.D. and CAMRA Consultants LLC.                |
# CAMRA Consultants LLC. weirmarkh@gmail.com, camraconsultants@gmail.com; weir.95@osu.edu                      |
#==============================================================================================================|

source("TriRand.r")
treatmentDist <- function(maxiter){
organisms <- c('crypto','giardia','rota','campy','eColi')
# Make the matrices -------------
coag <- matrix(NaN, nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))

filt_nc <- matrix(NaN, nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
filt_inline <- matrix(NaN, nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
filt_coag <- matrix(NaN, nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
filt_sand <- matrix(NaN, nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
filt_micro <- matrix(NaN, nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
filt_ultra <- matrix(NaN, nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))

bio_filt <- matrix(NaN, nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))

chlorine <- matrix(NaN,nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
chloromines <- matrix(NaN, nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
ozone <- matrix(NaN, nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))
chlorine_dioxide <- matrix(NaN,nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))

uv <- matrix(NaN, nrow=maxiter, ncol=length(organisms), dimnames = list(NULL, organisms))


# Simulate the reductions ---------------------
for(i in 1:maxiter)
{
  coag[,'crypto'][i] <- TriRand(0.46,1.8,3.78); 
  coag[,'giardia'][i] <- TriRand(0.48,1.3,0.48)
  coag[,'rota'][i] <- TriRand(0.34,1.7,4.06); 
  coag[,'campy'][i] <- TriRand(0.6,1.4,3.7); 
  coag[,'eColi'][i] <- TriRand(0.6,1.4,3.7)
  
  filt_nc[,'crypto'][i] <- TriRand(0.16,1.1,2.22); 
  filt_nc[,'giardia'][i] <- TriRand(0.07,1.0,2.39); 
  filt_nc[,'rota'][i] <- TriRand(0.11,0.6,1.59); 
  filt_nc[,'campy'][i] <- TriRand(0.2,0.4,1); 
  filt_nc[,'eColi'][i] <- TriRand(0.2,0.4,1)
  
  filt_inline[,'crypto'][i] <- TriRand(0.89,2.9,5.255); 
  filt_inline[,'giardia'][i] <- TriRand(1.3,3,4.6); 
  filt_inline[,'rota'][i] <- TriRand(0.16,0.3,1.5); 
  filt_inline[,'campy'][i] <- TriRand(0.86,1.4,1.95); 
  filt_inline[,'eColi'][i] <- TriRand(0.86,1.4,1.95)
  
  filt_coag[,'crypto'][i] <- TriRand(1.055,2.1,4.565); 
  filt_coag[,'giardia'][i] <- TriRand(0.72,1.7,3.75); 
  filt_coag[,'rota'][i] <- TriRand(0.265,0.8,2.11); 
  filt_coag[,'campy'][i] <- TriRand(0.44,0.9,1.42); 
  filt_coag[,'eColi'][i] <- TriRand(0.44,0.9,1.42)
  
  filt_sand[,'crypto'][i] <- TriRand(3.15,5,6.3); 
  filt_sand[,'giardia'][i] <- TriRand(4.07,4.7,5.87); 
  filt_sand[,'campy'][i] <- TriRand(0.66,2.1,3.7); 
  filt_sand[,'rota'][i] <- TriRand(1.288,2.4,4.64); 
  filt_sand[,'eColi'][i] <- TriRand(1.288,2.4,4.64)
  
  filt_micro[,'crypto'][i] <- TriRand(4.7,6.3,6.86); 
  filt_micro[,'giardia'][i] <- TriRand(5.95,6.8,6.975); 
  filt_micro[,'campy'][i] <- TriRand(0.1,0.5,2.4); 
  filt_micro[,'rota'][i] <- TriRand(2.86,3.9,6.66); 
  filt_micro[,'eColi'][i] <- TriRand(2.86,3.9,6.66)
  
  filt_ultra[,'crypto'][i] <- TriRand(5.74,6.4,6.96); 
  filt_ultra[,'giardia'][i] <- TriRand(5.04,6.55,6.965); 
  filt_ultra[,'rota'][i] <- TriRand(2.06,4.3,5.71); 
  filt_ultra[,'campy'][i] <- TriRand(8,8,8); 
  filt_ultra[,'eColi'][i] <- TriRand(8,8,8)
  
}

# Disinfection are point estimates ---------

chlorine[,'crypto'] <- runif(maxiter,0.02,0.02); 
chlorine[,'giardia'] <- runif(maxiter,1.65,1.65);
chlorine[,'rota'] <- runif(maxiter,8,8); 
chlorine[,'campy'] <- runif(maxiter,8,8) 
chlorine[,'eColi'] <- runif(maxiter,8,8)

chloromines[,'crypto'] <- runif(maxiter,0,0); 
chloromines[,'giardia'] <- runif(maxiter,0.09,0.09)
chloromines[,'rota'] <- runif(maxiter,0.7,0.7); 
chloromines[,'campy'] <- runif(maxiter,4.72,4.72)
chloromines[,'eColi'] <- runif(maxiter,1.36,1.36)

ozone[,'crypto'] <- runif(maxiter,3.01,3.01); 
ozone[,'giardia'] <- runif(maxiter,4,4)
ozone[,'rota'] <- runif(maxiter,4,4); 
ozone[,'campy'] <- runif(maxiter,8,8)
ozone[,'eColi'] <- runif(maxiter,8,8)

chlorine_dioxide[,'crypto'] <- runif(maxiter,0.2,0.2); 
chlorine_dioxide[,'giardia'] <- runif(maxiter,3.47,3.47);
chlorine_dioxide[,'rota'] <- runif(maxiter,8,8); 
chlorine_dioxide[,'campy'] <- runif(maxiter,4,4)
chlorine_dioxide[,'eColi'] <- runif(maxiter,4,4)

uv[,'crypto'] <- runif(maxiter, 0.966, 0.966)
uv[,'giardia'] <- runif(maxiter, 0.966, 0.966)
uv[,'rota'] <- runif(maxiter, 0.966, 0.966)
uv[,'campy'] <- runif(maxiter, 0.966, 0.966)
uv[,'eColi'] <- runif(maxiter, 0.966, 0.966)

# Biofilter is a point estimate
bio_filt[,'crypto'] <- runif(maxiter, 0.301, 0.301)
bio_filt[,'giardia'] <- runif(maxiter, 0.301, 0.301)
bio_filt[,'rota'] <- runif(maxiter, 0.301, 0.301)
bio_filt[,'campy'] <- runif(maxiter, 0.301, 0.301)
bio_filt[,'eColi'] <- runif(maxiter, 0.301, 0.301)

return(list(coag=coag,
            filt_nc=filt_nc,
            filt_coag=filt_coag,
            filt_inline=filt_inline,
            filt_sand=filt_sand,
            filt_micro=filt_micro,
            filt_ultra=filt_ultra,
            bio_filt=bio_filt,
            uv=uv,
            chlorine=chlorine,
            chloromines=chloromines,
            ozone=ozone,
            chlorine_dioxide=chlorine_dioxide
            )
       )
}
