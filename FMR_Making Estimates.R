
rm(list = ls(all = TRUE))

# setwd("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis")

####  Palette URL: http://paletton.com/#uid=72Y0u0khy-e7yT7cSJTlxtFqRqX
# Primary color:
pal.green  <- c("#B6EDCC", "#82D9A5", "#57C182", "#38AA66", "#199A4D")
# Secondary color (1):
pal.blue   <- c("#B4D8E8", "#7EB5CD", "#5392AD", "#357895", "#1B6787")
# Secondary color (2):
pal.orange <- c("#FFE3C3", "#FFCF98", "#FFBD73", "#ECA14D", "#D78223")
# Complement color:
pal.red    <- c("#FFCEC3", "#FFAC98", "#FF8E73", "#EC6C4D", "#D74523")

# install.packages(c("xlsx", "ape", "MCMCglmm"))

library(xlsx)
library(ape)
library(MCMCglmm)

# Sort data out:

data <- read.xlsx("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/Breeding FMR - Meta Analysis.xlsx", 1)

data$animal <- data$Match_name
data$Colony_explicit <- as.numeric(data$Colony_explicit)
data$log_Colony <- log10(data$Colony_explicit)

tr <- read.tree("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/bird tree seabird tree.phy")

sptree <- makeNodeLabel(tr, method = "number", prefix = "node")
treeAinv<-inverseA(sptree, nodes="TIPS")$Ainv

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.5.rda")

# plot(model2.5)
summary(model2.5)

################################################################################

## Make Estimates

# Example:
# White-faced storm petrel (Pelagodroma marina) 207 g

par(mar = c(4,4,2,0)+0.5)

# Incubation
posterior.mode(model2.5$Sol[,"(Intercept)"] + model2.5$Sol[,"animal.Pelagodroma_marina"] + mean(data$Lat)*model2.5$Sol[,"Lat"]+
                 log10(207)*model2.5$Sol[,"log_Mass"] + model2.5$Sol[,"PhaseIncubation"])
HPDinterval(model2.5$Sol[,"(Intercept)"] + model2.5$Sol[,"animal.Pelagodroma_marina"] + mean(data$Lat)*model2.5$Sol[,"Lat"]+
              log10(207)*model2.5$Sol[,"log_Mass"] + model2.5$Sol[,"PhaseIncubation"])

# Brood

posterior.mode(model2.5$Sol[,"(Intercept)"] + model2.5$Sol[,"animal.Pelagodroma_marina"] + mean(data$Lat)*model2.5$Sol[,"Lat"]+
                 log10(207)*model2.5$Sol[,"log_Mass"])
HPDinterval(model2.5$Sol[,"(Intercept)"] + model2.5$Sol[,"animal.Pelagodroma_marina"] + mean(data$Lat)*model2.5$Sol[,"Lat"]+
              log10(207)*model2.5$Sol[,"log_Mass"])

# Creche
posterior.mode(model2.5$Sol[,"(Intercept)"] + model2.5$Sol[,"animal.Pelagodroma_marina"] + mean(data$Lat)*model2.5$Sol[,"Lat"]+
                 log10(207)*model2.5$Sol[,"log_Mass"] + model2.5$Sol[,"PhaseCreche"])
HPDinterval(model2.5$Sol[,"(Intercept)"] + model2.5$Sol[,"animal.Pelagodroma_marina"] + mean(data$Lat)*model2.5$Sol[,"Lat"]+
              log10(207)*model2.5$Sol[,"log_Mass"] + model2.5$Sol[,"PhaseCreche"])


# Plot results:

par(mfrow = c(1,1), mar = c(4,4,0,0)+0.5)
plot(X.FMR ~ Mass_g, data = data, type = "n", axes = F, log = "xy", xlab = "Mass (g)", ylab = "FMR (kJ/ Day)")

axis(side = 1, las = 1, lwd = 1.5)
axis(side = 2, las = 1, lwd = 1.5)
box(lwd = 1.5)

lines(x = range(subset(data, Phase == "Incubation")$Mass_g),
      y = 10^(summary(model2.5)$solutions["(Intercept)","post.mean"]+
                mean(data$Lat)*summary(model2.5)$solutions["Lat","post.mean"]+
                log10(range(subset(data, Phase == "Incubation")$Mass_g))*
                summary(model2.5)$solutions["log_Mass","post.mean"]+
                summary(model2.5)$solutions["PhaseIncubation","post.mean"]),
      lwd = 1.5, col = pal.blue[5])
lines(x = range(subset(data, Phase == "Brood")$Mass_g),
      y = 10^(summary(model2.5)$solutions["(Intercept)","post.mean"]+
                mean(data$Lat)*summary(model2.5)$solutions["Lat","post.mean"]+
                log10(range(subset(data, Phase == "Brood")$Mass_g))*
                summary(model2.5)$solutions["log_Mass","post.mean"]),
      lwd = 1.5, col = pal.green[5], lty = 2)
lines(x = range(subset(data, Phase == "Creche")$Mass_g),
      y = 10^(summary(model2.5)$solutions["(Intercept)","post.mean"]+
                mean(data$Lat)*summary(model2.5)$solutions["Lat","post.mean"]+
                log10(range(subset(data, Phase == "Creche")$Mass_g))*
                summary(model2.5)$solutions["log_Mass","post.mean"]+
                summary(model2.5)$solutions["PhaseCreche","post.mean"]),
      lwd = 1.5, col = pal.red[5], lty = 3)


points(X.FMR ~ Mass_g, data = subset(data, Phase == "Incubation"),
       pch = 21, bg = "white", col = pal.blue[5], cex = 1)
points(X.FMR ~ Mass_g, data = subset(data, Phase == "Brood"),
       pch = 22, bg = "white", col = pal.green[5], cex = 1)
points(X.FMR ~ Mass_g, data = subset(data, Phase == "Creche"),
       pch = 23, bg = "white", col = pal.red[5], cex = 1)

# Add new points to plot:

points(x = 207,
       y = 10^(posterior.mode(model2.5$Sol[,"(Intercept)"] + model2.5$Sol[,"animal.Pelagodroma_marina"] + mean(data$Lat)*model2.5$Sol[,"Lat"]+
                                log10(207)*model2.5$Sol[,"log_Mass"] + model2.5$Sol[,"PhaseIncubation"])),
       pch = 22, col = "white", bg = pal.blue[5], cex = 1.5)

points(x = 207,
       y = 10^(posterior.mode(model2.5$Sol[,"(Intercept)"] + model2.5$Sol[,"animal.Pelagodroma_marina"] + mean(data$Lat)*model2.5$Sol[,"Lat"]+
                                log10(207)*model2.5$Sol[,"log_Mass"])),
       pch = 22, col = "white", bg = pal.green[5], cex = 1.5)


points(x = 207,
       y = 10^(posterior.mode(model2.5$Sol[,"(Intercept)"] + model2.5$Sol[,"animal.Pelagodroma_marina"] + mean(data$Lat)*model2.5$Sol[,"Lat"]+
                                log10(207)*model2.5$Sol[,"log_Mass"] + model2.5$Sol[,"PhaseCreche"])),
       pch = 22, col = "white", bg = pal.red[5], cex = 1.5)

legend(50, 6000, c("Creche", "Brood", "Incubation"),
       lty=c(1,2,3), # gives the legend appropriate symbols (lines)
       lwd=c(1.5,1.5,1.3),col=c(pal.red[5],pal.green[5],pal.blue[5])) # gives the legend lines the correct color and width

################################################################################

