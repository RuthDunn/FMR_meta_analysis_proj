rm(list = ls(all = TRUE))

library(rsconnect)
library(shiny)
library(ape)
library(MCMCglmm)
library(xlsx)
library(plyr)
library(shinydashboard)
library(stringr)

######################################################################

# Load in what we need outside the ui and server functions
# This means that they're only loaded once

# setwd(dir = "myapp/")

tr <- read.tree("data/bird tree seabird tree.phy")
sptree <- makeNodeLabel(tr, method = "number", prefix = "node")
seabirds <- read.csv("data/SeabirdSpecies.csv")
species_choices <- as.data.frame(sort(sptree$tip.label))
names(species_choices)[1] <- "animal"
species_choices <- merge(species_choices, seabirds, by = "animal")
species_choices<-species_choices[,c(1,2,4,6)]
rm(sptree, tr, seabirds)

species_choices_list <- split(species_choices$selector_caption, species_choices$Family_Common)

#

load(file = "data/model4.rda")

#####################################################################
#####################################################################

# Build User interface

# Input functions:

ui <- dashboardPage(
  
  skin = "black",
  
  dashboardHeader(title = "Seabird FMR Calculator",
                  titleWidth = 300),
  
  dashboardSidebar(disable = TRUE),
  
  dashboardBody(
    
          fluidRow(
            
            # Box 1 - Title, Authors & Info
            
            box(width = 4,
            background = "light-blue",
            status = "primary",
            
            h3("A model to estimate seabird field metabolic rates"),
            
            h5("Dunn, R.E., White, C.R., Green, J.A."),
            
            p("The",
            strong("Seabird FMR Calculator"),
            "is a web-based app which can be utilised to generate estimates of field metabolic rate (FMR)
              for any population of breeding seabird.",
            
            p("Daily FMR estimates are based on the outputs of a 
              phylogenetically informed meta-analytical model exploring the large-scale determinants
              of seabird FMR during the breeding season. The app requires inputs of species, bird mass,
              colony latitude and breeding phase. In return it generates an estimate of daily FMR alongside
              HPD confidence intervals.",
            
            p("We encourage the use of outputs generated from the",
            strong("Seabird FMR Calculator"),
                   "being utilised to complement future behavioural studies and to
            increase understanding of how energetic demands influence the role of seabirds
            as driving components of marine systems.
            ")))),
            
            # Box 2 - Selecting Model Inputs
            
            box(width = 4,
            status = "primary",
            h3("Model Inputs"),
          
            p("Enter model inputs below and then click",
              strong("Estimate FMR"),
              "to generate FMR estimates."),
            
            selectInput("species", label = "Select/ Type Species",
                      choices = species_choices_list),
            
            sliderInput("lat", label = "Colony Latitude (north/ south)",
                      min = 0, max = 90, value = 45),
            
            numericInput("mass", label = "Mean Bird Mass (g)",
                   value = 100),
            
            selectInput("phase", label = "Select Breeding Phase",
                  choices = list("Incubation",  #  Incubation is the intercept, so doesn't work
                                 "Brood",
                                 "Creche"),
                  selected = 1),
            
            actionButton("update" ,"Estimate FMR", icon("spinner"),
                   class = "btn btn-primary")),
            
            # Box 3 - Model Estimates
            
            box(width = 4,
                
                textOutput("your_model", h3),
                textOutput("selected_species"),
                textOutput("selected_lat"),
                textOutput("selected_mass"),
                textOutput("selected_phase"),
                
                hr(),
                
                textOutput("your_output", h3),
                textOutput("mode"),
                textOutput("lower_interval"),
                textOutput("upper_interval"))
)))

#####################################################################
#####################################################################

# Define server logic:

server <- function(input, output) {ntext_species <- selected_species <- eventReactive(input$update,{input$n})

ntext_lat <- selected_lat <- eventReactive(input$update, {input$n})

ntext_mass <- selected_mass <- eventReactive(input$update, {input$n})

ntext_phase <- selected_phase <- eventReactive(input$update, {input$n})

ntext_mode <- mode <- eventReactive(input$update, {input$n})

ntext_lower_interval <- lower_interval <- eventReactive(input$update, {input$n})
ntext_upper_interval <- upper_intervatl <- eventReactive(input$update, {input$n})

ntext_model <- your_model <- eventReactive(input$update, {input$n})
ntext_output <- your_output <- eventReactive(input$update, {input$n})

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

output$your_model <- renderText({
  ntext_model()
  isolate("Your Model:")
})

output$mode <-  renderText({
  ntext_mode()
  isolate(if(input$phase == "Brood"){
    paste("Estimate = ",
          round(10^(posterior.mode(model4$Sol[,"(Intercept)"] +
                                     model4$Sol[,paste("animal.", word(input$species,1), "_", word(input$species,2), sep = "")] +
                                     input$lat*model4$Sol[,"Lat"]+
                                     log10(input$mass)*model4$Sol[,"log_Mass"])), digits = 2),
          "kJ/Day")
  }else{
    paste("Estimate = ",
          round(10^(posterior.mode(model4$Sol[,"(Intercept)"] +
                                     model4$Sol[,paste("animal.", word(input$species,1), "_", word(input$species,2), sep = "")] +
                                     input$lat*model4$Sol[,"Lat"]+
                                     log10(input$mass)*model4$Sol[,"log_Mass"] +
                                     model4$Sol[,paste("Phase", input$phase, sep = "")])), digits = 2),
          "kJ/Day")
  }
  )})  

output$your_output <- renderText({
  ntext_output()
  isolate(paste("Your Output:"))
})

output$lower_interval <- renderText({
  ntext_lower_interval()
  isolate(if(input$phase == "Brood"){
    paste("Lower confidence interval:",
          round(10^(HPDinterval(model4$Sol[,"(Intercept)"] +
                                  model4$Sol[,paste("animal.", word(input$species,1), "_", word(input$species,2), sep = "")] +
                                  input$lat*model4$Sol[,"Lat"]+
                                  log10(input$mass)*model4$Sol[,"log_Mass"]))[1], digits = 2),
          "kJ/Day")
  }else{
    paste("Lower confidence interval:",
          round(10^(HPDinterval(model4$Sol[,"(Intercept)"] +
                                  model4$Sol[,paste("animal.", word(input$species,1), "_", word(input$species,2), sep = "")] +
                                  input$lat*model4$Sol[,"Lat"]+
                                  log10(input$mass)*model4$Sol[,"log_Mass"] +
                                  model4$Sol[,paste("Phase", input$phase, sep = "")]))[1], digits = 2),
          "kJ/ Day")
  }
  )})

output$upper_interval <- renderText({
  ntext_upper_interval()
  isolate(if(input$phase == "Brood"){
    paste("Upper confidence interval:",
          round(10^(HPDinterval(model4$Sol[,"(Intercept)"] +
                                  model4$Sol[,paste("animal.", word(input$species,1), "_", word(input$species,2), sep = "")] +
                                  input$lat*model4$Sol[,"Lat"]+
                                  log10(input$mass)*model4$Sol[,"log_Mass"]))[2], digits = 2),
          "kJ/Day")
  }else{
    paste("Upper confidence interval:",
          round(10^(HPDinterval(model4$Sol[,"(Intercept)"] +
                                  model4$Sol[,paste("animal.", word(input$species,1), "_", word(input$species,2), sep = "")] +
                                  input$lat*model4$Sol[,"Lat"]+
                                  log10(input$mass)*model4$Sol[,"log_Mass"] +
                                  model4$Sol[,paste("Phase", input$phase, sep = "")]))[2], digits = 2),
          "kJ/ Day")
  }
  )})}

#####################################################################
#####################################################################

# Run the app:

shinyApp(ui, server)
