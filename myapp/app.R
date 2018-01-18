rm(list = ls(all = TRUE))

library(rsconnect)
library(shiny)
library(ape)
library(MCMCglmm)
library(xlsx)

######################################################################

# Load in what we need outside the ui and server functions
# This means that they're only loaded once

# setwd(dir = "myapp/")

tr <- read.tree("data/bird tree seabird tree.phy")
sptree <- makeNodeLabel(tr, method = "number", prefix = "node")
species_choices <- as.list(sort(sptree$tip.label))
load(file = "data/model2.5c.rda")
plotcols <- c("#199A4D", "#1B6787", "#D78223")

#####################################################################

# Add elements to app:
# Input functions

ui <- fluidPage(

    titlePanel("Seabird FMR Calculator"),
                
    sidebarPanel(
      
      h4("Model Inputs"),
      
      p("Enter model inputs below and then click",
        strong("Update View"),
        "to calculate FMR estimates."),
      
      selectInput("species", label = "Select Species",
                  choices = species_choices),
      
      sliderInput("lat", label = "Colony Latitude (north/ south)",
                  min = 0, max = 90, value = 50),
      
      numericInput("mass", label = "Mean Bird Mass (g)",
                   value = 100),
      
      selectInput("phase", label = "Select Breeding Phase",
                  choices = list("Incubation",  #  Incubation is the intercept, so doesn't work
                                 "Brood",
                                 "Creche"),
                  selected = 1),
      
      actionButton("update" ,"Estimate FMR", icon("spinner"),
                   class = "btn btn-primary")
      
      ),
    
#####################################################################   

# Add elements to app:
# Output functions

    mainPanel(
      
      h4("A Model to Estimate Metabolic Rate in Seabirds"),
      
      h5("Dunn, R.E., White, C.R., Green, J.A."),
      
      p("The",
        strong("Seabird FMR Calculator"),
        "is the output of a series of phylogenetically-controlled meta-analytic mixed 
        effects models. These models incoorporate easily attainable 
        ecological and physiological parameters in order to make predictions 
        of FMR for seabird species where this has not previously been calculated."),
      
      h4("Your Model:"),
      
      textOutput("selected_species"),
      textOutput("selected_lat"),
      textOutput("selected_mass"),
      textOutput("selected_phase"),
      
      h4("Your Output:"),
      
      textOutput("mode"),
      textOutput("intervals"),
      
      plotOutput("plot")
      
    )
  )

#####################################################################

# Define server logic:

server <- function(input, output) {
  
  ntext_species <- selected_species <- eventReactive(input$update,{input$n})

  ntext_lat <- selected_lat <- eventReactive(input$update, {input$n})

  ntext_mass <- selected_mass <- eventReactive(input$update, {input$n})

  ntext_phase <- selected_phase <- eventReactive(input$update, {input$n})

  ntext_mode <- mode <- eventReactive(input$update, {input$n})

  ntext_intervals <- intervals <- eventReactive(input$update, {input$n})
  
    output$selected_species <- renderText({
    ntext_species()
    input$update
    isolate(paste("You have selected:", input$species))
  })
  
  output$selected_lat <- renderText({
    ntext_lat()
    input$update
    isolate(paste("You have identified a colony latitude of", input$lat, "N/ S"))
  })
  
  output$selected_mass <- renderText({
    ntext_mass()
    isolate(paste("You have chosen a mass of", input$mass, "g"))
  })
  
  output$selected_phase <- renderText({
    ntext_phase()
    isolate(paste("You have selected breeding phase:", input$phase))
  })
    
  output$mode <-  renderText({
    ntext_mode()
    isolate(if(input$phase == "Brood"){
      paste("Estimate = ",
          round(10^(posterior.mode(model2.5$Sol[,"(Intercept)"] +
                           model2.5$Sol[,paste("animal.", input$species, sep = "")] +
                           input$lat*model2.5$Sol[,"Lat"]+
                           log10(input$mass)*model2.5$Sol[,"log_Mass"])), digits = 2),
          "kJ/Day")
    }else{
      paste("Estimate = ",
            round(10^(posterior.mode(model2.5$Sol[,"(Intercept)"] +
                                 model2.5$Sol[,paste("animal.", input$species, sep = "")] +
                                 input$lat*model2.5$Sol[,"Lat"]+
                                 log10(input$mass)*model2.5$Sol[,"log_Mass"] +
                                 model2.5$Sol[,paste("Phase", input$phase, sep = "")])), digits = 2),
            "kJ/Day")
    }
  )})  
    
  output$intervals <- renderText({
    ntext_intervals()
    isolate(if(input$phase == "Brood"){
      paste("Confidence interval:",
          round(10^(HPDinterval(model2.5$Sol[,"(Intercept)"] +
                        model2.5$Sol[,paste("animal.", input$species, sep = "")] +
                        input$lat*model2.5$Sol[,"Lat"]+
                        log10(input$mass)*model2.5$Sol[,"log_Mass"])), digits = 2),
          "kJ/Day")
    }else{
      paste("Confidence interval:",
           round(10^(HPDinterval(model2.5$Sol[,"(Intercept)"] +
                              model2.5$Sol[,paste("animal.", input$species, sep = "")] +
                              input$lat*model2.5$Sol[,"Lat"]+
                              log10(input$mass)*model2.5$Sol[,"log_Mass"] +
                              model2.5$Sol[,paste("Phase", input$phase, sep = "")])), digits = 2),
           "kJ/ Day")
    }
    )})
}


# Run the app:

shinyApp(ui = ui, server = server)
