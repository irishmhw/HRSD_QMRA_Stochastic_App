# TriRand.r is an R (www.r-project.org) source code that is targeted at developing truncated normal distributed |
# random variables for use in R source codes.  The basics of the code are based on the workings of the normal   |
# distribution that is part of the R based package.                                                             |
# tnormal.r was written in entirety by Mark H. Weir EIT. Ph.D. with CAMRA Consultants LLC.    					|
#																												|
# msm package must be installed for the script to execute, 														|
# Execute the code and the function rtcauchy will be available, the syntax is as follows						|
# TriRand(minimum,likeliest,maximum)															  		        |
#																												|
# All use and reproduction rights are reserved by Mark H. Weir EIT. Ph.D. and CAMRA Consultants LLC.            |
# CAMRA Consultants LLC. 001-570-460-8459, weirmarkh@gmail.com, camraconsultants@gmail.com                      |
#===============================================================================================================|

require(stats)
TriRand <- function(minValue, likeValue, maxValue)
	{
			z = runif(1); .Random.seed[1:1];
			t = sqrt(z*(maxValue-minValue)*(likeValue-minValue))+minValue
			tt = maxValue-sqrt((1-z)*(maxValue-minValue)*(maxValue-likeValue))
			if (tt < likeValue) {return(t)} else{return(tt)}
	}

rtri <- function(n, minValue, likeValue, maxValue){
  return( sapply(rep(minValue, n), TriRand, likeValue, maxValue) )
}