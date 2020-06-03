# HRSD QMRA Stochastic Model and App

# App testing
Using a RStudio on a machine that has R version 3.6 or later, run the following snipets below. This will launch the app, and allow you to test it. You can then interact with the app in the screen that appears on a Mac, or in the browser on PC. You can also select to use the app in a browser in either platform. It is recommended to copy and paste the block of code all at one into RStudio. 

```r
# Copy and paste the entirety of this code chunk into RStudio console and press Enter

if("shiny" %in% rownames(installed.packages()) == FALSE) {install.packages("shiny", dependencies = TRUE); require(shiny)} else{require(shiny)}

if("shinyjs" %in% rownames(installed.packages()) == FALSE) {install.packages("shinyjs", dependencies = TRUE); require(shinyjs)} else{require(shinyjs)}

runGitHub('HRSD_QMRA_Stochastic_App', 'irishmhw')
```

# Purpose
  The app that is the primary reason for this repository. This app is meant for HRSD personnel to be able to interact with and use the stochastic quantitative microbial risk analysis (QMRA) model. The app is built using Shiny and therefore can be downlaoded in entirety from this repository and by using RStudio can be easily used. This will be the primary use of the model until OSU returns from COVID-19 shutdown to allow access to an app hosting server for this work. 
  
# Methods
As outlined in the final report, the QMRA model that the app is developed around is a stochastic QMRA model. This is the first of a kind QMRA application in two senses. First is the use of stochastic methods with the support to use your system's data, and also on a free software environment such as R. 

# Inputs
The user can enter exposure factors such as ingestion rate values in the spaces provided. The other inputs allow the the entry of a mean and standard deviation concentration of pathogen or monitored microorganism. However, more accuracy will be provided by uploading data as a .csv file in the provided upload dialog boxes. After upoloading or entering the data desired, then the efficiency values of the processes can be set using the sliders in the app, setting the number of iterations and clicking run. 

For testing purposes it is recommended to run 1000 iterations due to overall app run time for testing and evaluation. At 1,000 iterations the run time should be ~20 seconds depending on local computer resources, at 5,000 iterations the run time should be ~3 minutes. There is no appreciable difference in outcomes between these two iteration options, therefore, 1,000 iterations is recommended for testing purposes. 

# Output
Violin plots and boxplots are visualization outputs that are seen to the right of the entry sidebar. These are representations of the density or interquartile range (violin and boxplot respectively) of the risk estiamtes


# Updates
