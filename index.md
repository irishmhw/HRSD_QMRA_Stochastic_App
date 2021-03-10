# Stochastic QMRA Model Application for HRSD 

## Purpose
This app was developed in the shiny interface in R, but outside of R-studio, for more information on that please see the use section below. The app was developed under a contract with the Hampton Roads Sanitation District (HRSD) for their SWIFT wastewater reuse facility in VA, USA. This facility is a state-of-the-art wastewater reuse system that was also built to be an educational setting for Environmental Engineering and Science. In an effort to support their decision making and risk communication with the communities downfield of their discharge into the groundwater, a quantitative microbial risk assessment (QMRA) model was developed. Contracting with Mark H. Weir Ph.D. ([irishmhw](https://github.com/irishmhw)) and [Nikki Beestch](https://www.linkedin.com/in/nikki-beetsch-ahrns-247ab89/) with NSF International, this QMRA model was developed. To facilitate the use of the model the shiny app that can be found at the repository link or can be downlodaed too. 

## Use of the QMRA Model
Follow the instructions in the list below and you will be able to operate the app. You __will require__ R (tested last on version 3.4.3). The more updated the version of R the more likely that dependencies for shiny apps will not be present. Downgrade your R implementation or conduct a fresh install on a previous stable version of R on another or virtual machine if needed.

- Download the repository using the green code button on the repo.
- unzip the file
- In your R console of choosing enter the follwoing without the bullet
    - ``` if("shiny" %in% row.names(installed.packages())==FALSE){install.packages("shiny"); require(shiny)}else{require(shiny)} ``` and ```if("shinyjs" %in% rownames(installed.packages()) == FALSE) {install.packages("shinyjs"); require(shinyjs)}else{require(shinyjs)}```
    - This will automatically install the package and load it if you do not have it installed. This can be edited to inlcude the (..., depdendencies=TRUE) in the ```install.pacakges()``` statement, but that is not necesaarily neeed to operate the app. 
    - After the shiny package is installed and you have seen that it is loaded then
    - ```runGitHub('HRSD_QMRA_Stochastic_App', 'irishmhw')``` OR ```runApp("app.r")```
    - This will open a new browser window and you can interact with the QMRA app
- Questions can be delivered to Dr. Weir, please be patient for a response.

## Live Operation of the Model

Coming soon an implementation of the model here. 








___________ General GitHub Pages Information that I kept in ___________


For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/irishmhw/HRSD_QMRA_Stochastic_App/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
