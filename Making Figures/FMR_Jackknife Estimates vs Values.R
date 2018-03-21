
rm(list = ls(all = TRUE))

# install.packages(c("xlsx", "ape", "MCMCglmm", "gplots"))

library(xlsx)
library(ape)
library(MCMCglmm)

################################################################################

# setwd("G:/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Jackknifing/")
# file.names <- list.files()
# file.names.df <- as.data.frame(list.files())

data <- read.xlsx("M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Data/Breeding FMR - Meta Analysis.xlsx", 1)

data$animal <- data$Match_name
data$Colony_explicit <- as.numeric(data$Colony_explicit)
data$log_Colony <- log10(data$Colony_explicit)

################################################################################

All.estimates <- NULL

for (i in 1:nrow(file.names.df)){
  
  # i = 6
  
  print(i)
  
  model <- file.names[i]
  row <- data[i,]
  
  load(paste(file = "G:/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Jackknifing/",
             model, sep = ""))
  
  load(paste(file = "G:/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Jackknifing/mod_14.rda"))
  
  #
  
  if(row$Phase == "Incubation"){
  
  estimate <- posterior.mode(mod$Sol[,"(Intercept)"] +
                               mod$Sol[,paste("animal.",row$Match_name, sep = "")] +
                               mean(row$Lat)*mod$Sol[,"Lat"] +
                               log10(row$Mass_g)*mod$Sol[,"log_Mass"] +
                               mod$Sol[,paste("Phase", row$Phase, sep = "")])
  
  intervals <- HPDinterval(mod$Sol[,"(Intercept)"] +
                             mod$Sol[,paste("animal.",row$Match_name, sep = "")] +
                             mean(data$Lat)*mod$Sol[,"Lat"]+
                             log10(row$Mass_g)*mod$Sol[,"log_Mass"] +
                             mod$Sol[,paste("Phase",row$Phase, sep = "")])
  }
  
  if(row$Phase == "Brood"){
    
    estimate <- posterior.mode(mod$Sol[,"(Intercept)"] +
                                 mod$Sol[,paste("animal.",row$Match_name, sep = "")] +
                                 mean(row$Lat)*mod$Sol[,"Lat"] +
                                 log10(row$Mass_g)*mod$Sol[,"log_Mass"])
    
    intervals <- HPDinterval(mod$Sol[,"(Intercept)"] +
                               mod$Sol[,paste("animal.",row$Match_name, sep = "")] +
                               mean(data$Lat)*mod$Sol[,"Lat"] +
                               log10(row$Mass_g)*mod$Sol[,"log_Mass"])
    }
  
  if(row$Phase == "Creche"){
    
    estimate <- posterior.mode(mod$Sol[,"(Intercept)"] +
                                 mod$Sol[,paste("animal.",row$Match_name, sep = "")] +
                                 mean(row$Lat)*mod$Sol[,"Lat"] +
                                 log10(row$Mass_g)*mod$Sol[,"log_Mass"] +
                                 mod$Sol[,paste("Phase", row$Phase, sep = "")])
    
    intervals <- HPDinterval(mod$Sol[,"(Intercept)"] +
                               mod$Sol[,paste("animal.",row$Match_name, sep = "")] +
                               mean(data$Lat)*mod$Sol[,"Lat"]+
                               log10(row$Mass_g)*mod$Sol[,"log_Mass"] +
                               mod$Sol[,paste("Phase",row$Phase, sep = "")])
    }
  
  #
  
  All.estimates <- rbind(All.estimates, data.frame(model, estimate, intervals))
 
  # write.csv(All.estimates,
  #           "M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Jackknifing/Jackknife_Estimates.csv")
  
}

#




################################################################################

library(ggplot2)

All.estimates <- read.csv("M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Jackknifing/Jackknife_Estimates.csv")
data <- cbind(All.estimates, data)

# plot(Inc.estimates$estimate, Incubation$log_FMR)

# plotCI(x = data[1:13,]$log_FMR, y = All.estimates$estimate,
#        ui = All.estimates$upper,
#        li = All.estimates$lower,
#        xlab = "log FMR Extracted Value",
#        ylab = "log FMR Estimates", gap = FALSE)
# abline(0,1, lty = 2)

# Colour blind friendly palette:

hello <- theme_bw() %+replace% theme(panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank())

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot <- ggplot(data, aes(x = log_FMR, y = estimate)) +
  geom_point(aes(x = log_FMR, y = estimate, color = Phase)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.01, color = Phase)) +
  scale_colour_manual(values=cbPalette[c(2, 4, 6)]) +
  # scale_shape_manual(values = c(0:4, 6, 8, 15:17)) + 
  xlab("True log FMR Values (kJ/ day)")+
  ylab("Jackknife Estimate log FMR Values (kJ/ day)") +
  geom_abline(linetype = 2)

theme_set(hello)

plot

# ggsave("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Plots/jackknife_estimates.png",
#        width = 8, height = 5)
