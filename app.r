list.of.packages <- c('jsonlite','dplyr','reshape','ggplot2','gridExtra')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

require(shiny)
require(shinyjs)
require(jsonlite)
require(dplyr)


orgJson <- fromJSON('./organisms.json')
organisms <-
  orgJson$organisms[c('name',
                      'short',
                      'beta',
                      'optimizedParam',
                      'model',
                      'illGivenInfect',
                      'DALYs')]

# Use the concentration values supplied.
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

coag <- filt <- data.frame(
  name = c('coag'),
  title = c('Coag')
)
filt <- data.frame(
  name = c('filt_inline','filt_nc','filt_coag','filt_sand','filt_micro','filt_ultra'),
  title = c('Inline', 'No Coag', 'Coag', 'Sand', 'Micro', 'Ultra')
)

disinfect<- data.frame(
  name = c('chlorine','chloromines','ozone','chlorine_dioxide','uv'),
  title = c('Chlorine','Chloromines','Ozone','Chlorine Dioxide','UV')
)

ui <- fluidPage(
  # App title ----
  titlePanel("HRSD SWIFT QMRA"),
  sidebarPanel(
    actionButton(
      'run',
      "Run",
      icon = NULL,
      width = '100%',
      class = 'btn-success'
    ),
    div(
      h3('Treatment'), 
      selectInput('coag','Coagulation:', coag$title),
      sliderInput('coagEff', 'Efficiency', 0, 1, 1),
      selectInput('disinfect1', 'Disinfection 1: ', disinfect$title, selected='Ozone'),
      sliderInput('disinfect1Eff', 'Efficiency', 0, 1, 1),
      selectInput('bioFilt','Bio Filtration:', c('biofiltration')),
      sliderInput('bioFiltEff', 'Efficiency', 0, 1, 1),
      selectInput('filt','Filtration:', filt$title),
      sliderInput('filtEff', 'Efficiency', 0, 1, 1),
      selectInput('disinfect2', 'Disinfection 2: ', disinfect$title, selected='UV'),
      sliderInput('disinfect2Eff', 'Efficiency', 0, 1, 1),
      selectInput('disinfect3', 'Disinfection 3: ', disinfect$title, selected='Chlorine'),
      sliderInput('disinfect3Eff', 'Efficiency', 0, 1, 1)
    ),
    div(
      h3('Ingestion'),
      numericInput('ingest_min', 'Min:', 0.5, min = 0.001, max = 10),
      numericInput('ingest_mean', 'Likeliest:', 1, min = 0.001, max = 10),
      numericInput('ingest_max', 'Max:', 2.5, min = 0.001, max = 10)
    ),
    div(
      h3('Concentration'),
      numericInput('conc_mean', 'Mean:', 0.25, min = 0.001, max = 0.999, step=0.01),
      numericInput('conc_sd', 'Std Dev:', 0.22, min = 0.001, max = 0.999, step=0.01)
    ),
    fileInput("crypto_file", "Upload Crypto",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")
    ),
    fileInput("giardia_file", "Upload Giardia",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")
    ),
    fileInput("rota_file", "Upload Rotavirus",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")
    ),
    fileInput("campy_file", "Upload Campylobacter",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")
    ),
    fileInput("ecoli_file", "Upload Ecoli",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")
    )
  ),
  mainPanel(
    numericInput(
      'maxiter',
      label = 'Max Iterations',
      1000,
      min = 1000,
      max = 100000,
      step = 100
    ),
    h4(textOutput(outputId = 'serverMessage'), height=100),
    
    # Output: Histograms ----
    htmlOutput(outputId = "mcManifest"),
    uiOutput(outputId = "dataInput"),
    plotOutput(outputId = "treatPlot", height=3000),
    htmlOutput(outputId = "mcData"),
    plotOutput(outputId = "mcPlot", height=1500)
    
  )
)

server <- function(input, output, session) {
  source('monte_carlo.r')
  source('daly.r')
  
  observeEvent(input$crypto_file, {
    cryptoData <<- read.csv(input$crypto_file$datapath, header=TRUE, col.names=c('x'))
    output$dataInput <- renderUI(div(
      h2("Input Data"),
      renderDataTable(getDataInput())
    ))
  })
  
  observeEvent(input$giardia_file, {
    giardiaData <<- read.csv(input$giardia_file$datapath, header=TRUE, col.names=c('x'))
    output$dataInput <- renderUI(div(
      h2("Input Data"),
      renderDataTable(getDataInput())
    ))
  })
  
  observeEvent(input$rota_file, {
    rotaData <<- read.csv(input$rota_file$datapath, header=TRUE, col.names=c('x'))
    output$dataInput <- renderUI(div(
      h2("Input Data"),
      renderDataTable(getDataInput())
    ))
  })
      
  observeEvent(input$campy_file, {
    campyData <<- read.csv(input$campy_file$datapath, header=TRUE, col.names=c('x'))
    output$dataInput <- renderUI(div(
      h2("Input Data"),
      renderDataTable(getDataInput())
    ))
  })
  
  observeEvent(input$ecoli_file, {
    eColiData <<- read.csv(input$ecoli_file$datapath, header=TRUE, col.names=c('x'))
    output$dataInput <- renderUI(div(
      h2("Input Data"),
      renderDataTable(getDataInput())
    ))
  })
  
  observeEvent(input$run, {
      # Perform the Monte Carlo when the "Run" button is pressed.
      if(input$run){
        # Grab the other inputs. These provide the
        # additional arguments to the MC.
        inputData <- getDataInput() # Any user file data
        coagInp <- coag %>% filter(coag$title == input$coag)
        filtInp <- filt %>% filter(filt$title == input$filt)
        disinfect1Inp <- disinfect %>% filter(disinfect$title == input$disinfect1)
        disinfect2Inp <- disinfect %>% filter(disinfect$title == input$disinfect2)
        disinfect3Inp <- disinfect %>% filter(disinfect$title == input$disinfect3)
        
        efficiency <- list(
                           coag=input$coagEff,
                           disinfect1=input$disinfect1Eff,
                           filt=input$filtEff,
                           bioFilt=input$bioFiltEff,
                           disinfect2=input$disinfect2Eff,
                           disinfect3=input$disinfect3Eff
                           )
        
        # qmra.MonteCarlo returns a list of data.frames.
        # These are used for plotting and displaying data.
        mcRes <- qmra.MonteCarlo(inputData, 
                                 input$maxiter,
                                 coagName = coagInp$name,
                                 filtName = filtInp$name,
                                 disinfectName1 = disinfect1Inp$name,
                                 disinfectName2 = disinfect2Inp$name,
                                 efficiency = efficiency,
                                 ingest = list(min=input$ingest_min, mean=input$ingest_mean, max=input$ingest_max)
                                )
        
        annRisk <- data.frame(crypto = numeric(input$maxiter),
                             giardia = numeric(input$maxiter),
                             rota = numeric(input$maxiter),
                             campy = numeric(input$maxiter),
                             eColi = numeric(input$maxiter))
        
        dalys <- data.frame(crypto = numeric(input$maxiter),
                            giardia = numeric(input$maxiter),
                            rota = numeric(input$maxiter),
                            campy = numeric(input$maxiter),
                            ecoli = numeric(input$maxiter))
        
        for (i in 1:nrow(organisms)){
          org <- organisms[i,]
          # FIXME Do we need theses?
          annRisk[,org$short] <- AnnualRiskOfIllness(mcRes$riskMorb[, org$short], org$illGivenInfect)
          dalys[,org$short] <- calcDALYs(annRisk[,org$short], org$DALYs)
        }
  
      }
      
      # Handle the results.
      
      # First the data itself.
      prbs <- c(0.005, 0.05, 0.5, 0.95, 0.995)
      output$mcData <- renderUI(
        
        if(input$run){
          div(
            h1("Results"),
            h2("Risk"),
            renderTable(apply(mcRes$risk, 2, quantile, probs=prbs), rownames = TRUE, digits = -4),
            h2("Risk Morb"),
            renderTable(apply(mcRes$riskMorb, 2, quantile, probs=prbs), rownames = TRUE, digits = -4),
            h2("Annual Risk Morb"),
            renderTable(apply(annRisk, 2, quantile, probs=prbs), rownames = TRUE, digits = -4),
            h2("DALYs"),
            renderTable(apply(dalys, 2, quantile, probs=prbs), rownames = TRUE, digits = -4)
          )
        }
      )
      
      # Create the main plots (violin and box) for each organism
      output$treatPlot <- renderPlot(
        if(input$run){
          plots <- list()
          for (i in 1:nrow(organisms)){
            org <- organisms[i,]$short
            
            cols <- c(
              input$coag, 
              paste(input$disinfect1,'1',sep=''), 
              input$bioFilt, input$filt,  
              paste(input$disinfect2,'2',sep=''),  
              paste(input$disinfect3,'3',sep='')
            )
            dose <- data.frame(cbind(mcRes$postCoagDose[,org], mcRes$postDisinfect1Dose[,org],  mcRes$postBioFiltDose[,org],  mcRes$postFiltDose[,org],  mcRes$postDisinfect2Dose[,org], mcRes$postDisinfect3Dose[,org]))
            risk <- data.frame(cbind(mcRes$postCoagRisk[,org], mcRes$postDisinfect1Risk[,org], mcRes$postBioFiltRisk[,org], mcRes$postFiltRisk[,org],  mcRes$postDisinfect2Risk[,org], mcRes$postDisinfect3Risk[,org]))
            illness <- data.frame(cbind(mcRes$postCoagMorb[,org], mcRes$postDisinfect1Morb[,org], mcRes$postBioFiltMorb[,org], mcRes$postFiltMorb[,org],  mcRes$postDisinfect2Morb[,org], mcRes$postDisinfect3Morb[,org]))
            
            colnames(dose) <- cols
            colnames(risk) <- cols
            colnames(illness) <- cols
            
            plots[[i]] <- suppressWarnings(qmra.treatmentPlots(organisms[i,]$name, dose, risk, illness))
          }
          
          grid.arrange(grobs=plots, nrows=5, ncol=1, heights=unit(rep(8,5), rep('in', 5)))
        }
      )
      
      # These are additional plots used during development.
      # output$mcPlot <- renderPlot(
      #  if(input$run){
      #      qmra.MonteCarloPlot(mcRes)
      #  }
      # )
  })
}

shinyApp(ui, server)