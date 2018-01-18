
rm(list = ls(all = TRUE))

# install.packages(c("xlsx", "ape", "MCMCglmm"))

library(xlsx)
library(ape)
library(MCMCglmm)
library(gplots)

################################################################################

data <- read.xlsx("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/Breeding FMR - Meta Analysis.xlsx", 1)

data$animal <- data$Match_name
data$Colony_explicit <- as.numeric(data$Colony_explicit)
data$log_Colony <- log10(data$Colony_explicit)

# Split into Incubation, Brood and Creche

Incubation <- split(data, data$Phase)[[3]]
Incubation$Phase <- factor(Incubation$Phase)

Brood <- split(data, data$Phase)[[1]]
Brood$Phase <- factor(Brood$Phase)

Creche <- split(data, data$Phase)[[2]]
Creche$Phase <- factor(Creche$Phase)

################################################################################

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.5d.rda")

################################################################################

# Make estimates for every species and mass in database

par(mfrow = c(1,3))

# Incubation

Inc.estimates <- NULL

for (i in 1:nrow(Incubation)){
  print(i)

  # i = 13
  
estimate <- posterior.mode(model2.5$Sol[,"(Intercept)"] +
                             model2.5$Sol[,paste("animal.",Incubation$Match_name[i], sep = "")] +
                             mean(Incubation$Lat[i])*model2.5$Sol[,"Lat"] +
                             log10(Incubation$Mass_g[i])*model2.5$Sol[,"log_Mass"] +
                             model2.5$Sol[,paste("Phase","Incubation", sep = "")])

intervals <- HPDinterval(model2.5$Sol[,"(Intercept)"] +
                           model2.5$Sol[,paste("animal.",Incubation$Match_name[i], sep = "")] +
                           mean(data$Lat)*model2.5$Sol[,"Lat"]+
                           log10(Incubation$Mass_g[i])*model2.5$Sol[,"log_Mass"] +
                           model2.5$Sol[,paste("Phase","Incubation", sep = "")])

# print(estimate)

Inc.estimates <- rbind(Inc.estimates, data.frame(i, estimate, intervals))

}

##

# plot(Inc.estimates$estimate, Incubation$log_FMR)
plotCI(x = Incubation$log_FMR, y = Inc.estimates$estimate,
       ui = Inc.estimates$upper,
       li = Inc.estimates$lower,
       xlab = "Incubation log FMR Extracted Value",
       ylab = "Incubation log FMR Estimates", , gap = FALSE)
abline(0,1, lty = 2)

##

################################################################################

# Brood

Brood.estimates <- NULL

for (i in 1:nrow(Brood)){
  print(i)
  
  # i = 43
  
  estimate <- posterior.mode(model2.5$Sol[,"(Intercept)"] +
                               model2.5$Sol[,paste("animal.",Brood$Match_name[i], sep = "")] +
                               mean(Brood$Lat[i])*model2.5$Sol[,"Lat"] +
                               log10(Brood$Mass_g[i])*model2.5$Sol[,"log_Mass"])
  
  intervals <- HPDinterval(model2.5$Sol[,"(Intercept)"] +
                             model2.5$Sol[,paste("animal.",Brood$Match_name[i], sep = "")] +
                             mean(data$Lat)*model2.5$Sol[,"Lat"]+
                             log10(Brood$Mass_g[i])*model2.5$Sol[,"log_Mass"])
  
  # print(estimate)
  
  Brood.estimates <- rbind(Brood.estimates, data.frame(i, estimate, intervals))
  
}

##

# plot(Brood.estimates$estimate, Brood$log_FMR)
plotCI(x = Brood$log_FMR, y = Brood.estimates$estimate,
       ui = Brood.estimates$upper,
       li = Brood.estimates$lower,
       xlab = "Brood log FMR Estimates",
       ylab = "Brood log FMR Extracted Value", gap = FALSE)
abline(0,1, lty = 2)

##

################################################################################

# Creche

Creche.estimates <- NULL

for (i in 1:nrow(Creche)){
  print(i)
  
  # i = 43
  
  estimate <- posterior.mode(model2.5$Sol[,"(Intercept)"] +
                               model2.5$Sol[,paste("animal.",Creche$Match_name[i], sep = "")] +
                               mean(Creche$Lat[i])*model2.5$Sol[,"Lat"] +
                               log10(Creche$Mass_g[i])*model2.5$Sol[,"log_Mass"] +
                               model2.5$Sol[,paste("Phase","Creche", sep = "")])
  
  intervals <- HPDinterval(model2.5$Sol[,"(Intercept)"] +
                             model2.5$Sol[,paste("animal.",Creche$Match_name[i], sep = "")] +
                             mean(data$Lat)*model2.5$Sol[,"Lat"]+
                             log10(Creche$Mass_g[i])*model2.5$Sol[,"log_Mass"] +
                             model2.5$Sol[,paste("Phase","Creche", sep = "")])
  
  # print(estimate)
  
  Creche.estimates <- rbind(Creche.estimates, data.frame(i, estimate, intervals))
  
}

##

# plot(Creche.estimates$estimate, Creche$log_FMR)
plotCI(x = Creche$log_FMR, y = Creche.estimates$estimate,
       ui = Creche.estimates$upper,
       li = Creche.estimates$lower,
       xlab = "Creche log FMR Estimates",
       ylab = "Creche log FMR  Extracted Value", gap = FALSE)

abline(0,1, lty = 2)

##

################################################################################







################################################################################







################################################################################

# Without the log thing

par(mfrow = c(1,3))

# Incubation

Inc.estimates <- NULL

for (i in 1:nrow(Incubation)){
  print(i)
  
  # i = 13
  
  estimate <- 10^posterior.mode(model2.5$Sol[,"(Intercept)"] +
                               model2.5$Sol[,paste("animal.",Incubation$Match_name[i], sep = "")] +
                               mean(Incubation$Lat[i])*model2.5$Sol[,"Lat"] +
                               log10(Incubation$Mass_g[i])*model2.5$Sol[,"log_Mass"] +
                               model2.5$Sol[,paste("Phase","Incubation", sep = "")])
  
  intervals <- 10^HPDinterval(model2.5$Sol[,"(Intercept)"] +
                             model2.5$Sol[,paste("animal.",Incubation$Match_name[i], sep = "")] +
                             mean(data$Lat)*model2.5$Sol[,"Lat"]+
                             log10(Incubation$Mass_g[i])*model2.5$Sol[,"log_Mass"] +
                             model2.5$Sol[,paste("Phase","Incubation", sep = "")])
  
  # print(estimate)
  
  Inc.estimates <- rbind(Inc.estimates, data.frame(i, estimate, intervals))
  
}

##

# plot(Inc.estimates$estimate, Incubation$log_FMR)
plotCI(x = Incubation$X.FMR, y = Inc.estimates$estimate,
       ui = Inc.estimates$upper,
       li = Inc.estimates$lower,
       xlab = "Incubation FMR Extracted Value",
       ylab = "Incubation FMR Estimates",
       xlim = c(0, 7000),
       ylim = c(0, 10000))
abline(0,1, lty = 2)

##

################################################################################

# Brood

Brood.estimates <- NULL

for (i in 1:nrow(Brood)){
  print(i)
  
  # i = 43
  
  estimate <- 10^posterior.mode(model2.5$Sol[,"(Intercept)"] +
                               model2.5$Sol[,paste("animal.",Brood$Match_name[i], sep = "")] +
                               mean(Brood$Lat[i])*model2.5$Sol[,"Lat"] +
                               log10(Brood$Mass_g[i])*model2.5$Sol[,"log_Mass"])
  
  intervals <- 10^HPDinterval(model2.5$Sol[,"(Intercept)"] +
                             model2.5$Sol[,paste("animal.",Brood$Match_name[i], sep = "")] +
                             mean(data$Lat)*model2.5$Sol[,"Lat"]+
                             log10(Brood$Mass_g[i])*model2.5$Sol[,"log_Mass"])
  
  # print(estimate)
  
  Brood.estimates <- rbind(Brood.estimates, data.frame(i, estimate, intervals))
  
}

##

# plot(Brood.estimates$estimate, Brood$log_FMR)
plotCI(x = Brood$X.FMR, y = Brood.estimates$estimate,
       ui = Brood.estimates$upper,
       li = Brood.estimates$lower,
       xlab = "Brood FMR Extracted Value",
       ylab = "Brood FMR Estimates",
       xlim = c(0, 7000),
       ylim = c(0, 10000))
abline(0,1, lty = 2)

##

################################################################################

# Creche

Creche.estimates <- NULL

for (i in 1:nrow(Creche)){
  print(i)
  
  # i = 43
  
  estimate <- 10^posterior.mode(model2.5$Sol[,"(Intercept)"] +
                               model2.5$Sol[,paste("animal.",Creche$Match_name[i], sep = "")] +
                               mean(Creche$Lat[i])*model2.5$Sol[,"Lat"] +
                               log10(Creche$Mass_g[i])*model2.5$Sol[,"log_Mass"] +
                               model2.5$Sol[,paste("Phase","Creche", sep = "")])
  
  intervals <- 10^HPDinterval(model2.5$Sol[,"(Intercept)"] +
                             model2.5$Sol[,paste("animal.",Creche$Match_name[i], sep = "")] +
                             mean(data$Lat)*model2.5$Sol[,"Lat"]+
                             log10(Creche$Mass_g[i])*model2.5$Sol[,"log_Mass"] +
                             model2.5$Sol[,paste("Phase","Creche", sep = "")])
  
  # print(estimate)
  
  Creche.estimates <- rbind(Creche.estimates, data.frame(i, estimate, intervals))
  
}

##

# plot(Creche.estimates$estimate, Creche$log_FMR)
plotCI(x = Creche$X.FMR, y = Creche.estimates$estimate,
       ui = Creche.estimates$upper,
       li = Creche.estimates$lower,
       xlab = "Creche FMR Extracted Value",
       ylab = "Creche FMR Estimates",
       xlim = c(0, 7000),
       ylim = c(0, 10000))
abline(0,1, lty = 2)

##
