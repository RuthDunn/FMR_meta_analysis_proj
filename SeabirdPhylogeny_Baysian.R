
# rm(list = ls(all = TRUE))

#setwd("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis")
#setwd("~/Dropbox (White Lab)/People/Jon")

# install.packages("rncl")


####  Palette URL: http://paletton.com/#uid=72Y0u0khy-e7yT7cSJTlxtFqRqX
# Primary color:
pal.green  <- c("#B6EDCC", "#82D9A5", "#57C182", "#38AA66", "#199A4D")
# Secondary color (1):
pal.blue   <- c("#B4D8E8", "#7EB5CD", "#5392AD", "#357895", "#1B6787")
# Secondary color (2):
pal.orange <- c("#FFE3C3", "#FFCF98", "#FFBD73", "#ECA14D", "#D78223")
# Complement color:
pal.red    <- c("#FFCEC3", "#FFAC98", "#FF8E73", "#EC6C4D", "#D74523")


data <- read.csv("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Breeding_FMR_Physiology.csv")

data$animal <- as.character(data$animal)
data$Latitude <- abs(data$Latitude)
data$log_Colony <- log10(data$Colony)

library(rotl)
library(nlme)
library(lme4)
library(MCMCglmm)
library(ape)
library(car)
library(RODBC)
library(rncl)





      # model_simple <- MCMCglmm(phen ~ cofactor,
          # random = ~phylo,
          # family = "gaussian", ginverse = list(phylo = inv.phylo$Ainv), prior=prior,
          # data = data, nitt = 5000000, burnin = 1000, thin = 500)

# These are "flat", "noninformative" priors:
prior <- list(G = list(G1 = list(V = 1, nu = 0.001)),    # Random effect
            R = list(V = 1, nu = 0.001))                 # Fixed effect

# Model without phylogeny effect
model2.1<- MCMCglmm( log_FMR ~ log_FMR_Mass + Phase +  Latitude + log_Colony + Average_Brood,
                     random =~animal,
                     data = data,  prior = prior, 
                     burnin = 10000, nitt = 100000, thin = 200, verbose = T)

# Nitt = Total number of runs > 10,000
# BurnIn (running time before) = Recommended to be 10% of iterations and > 1000
# Thin = Sampling frequency

summary(model2.1)
plot(model2.1)
# Left hand plots should show no autocorrelation
# Right hand plots should show normal distributions


# Linear model of FMR and variables
m2.1 <- lm(log_FMR ~ log_FMR_Mass + Phase +  Latitude + log_Colony + Average_Brood, 
           data = data, na.action = na.omit)

# Plots of FMR against each variables, taking into account the effect of the other variables
avPlots(m2.1)




###
### Select species to be included within phylogeny analysis
###

# Red-legged kittiwake (Rissa brevirostris)
# Australian pelican (Pelecanus conspicillatus)
# Peruvian diving petrel (Pelecanoides garnotii)
# extra.species <- c("Rissa brevirostris", "Pelecanus conspicillatus", "Pelecanoides garnotii")
# resolved_names <- tnrs_match_names(c(data$animal,extra.species))

# Set this to TRUE to download the tree from ROTL
download.tree <- TRUE

# Use seabirds ott_ids (Open Tree Taxonomy identifier) to match and extract seabirds from the Open Tree tree
seabirds <- read.csv("SeabirdSpecies.csv")
seabirds$animal <- as.character(seabirds$animal)
# str(data)
resolved_names1 <- tnrs_match_names(head(seabirds$animal,250))
resolved_names2 <- tnrs_match_names(tail(seabirds$animal,101))
resolved_names3 <- tnrs_match_names(data$animal)
resolved_names <- rbind(resolved_names1, resolved_names2,resolved_names3)

# resolved_names <- tnrs_match_names(c(data$animal,seabirds$animal))
# Get the tree with just those tips:
tr <- tol_induced_subtree(ott_ids=ott_id(resolved_names))

# Open Tree trees have no branch lengths but MCMCglmm accepts topology as input

taxon_map <- structure (resolved_names$search_string, names = resolved_names$unique_name)

# Tree contains ID numbers which don't therefore match our data
# Use function to remove extra info at label tips
# Then use taxon map to replace the tip labels with sp. names from the data
otl_tips <- strip_ott_ids(tr$tip.label, remove_underscores=TRUE)
tr$tip.label <- taxon_map[ otl_tips ]
# Remove node labels by setting node.label to null
tr$node.label <- NULL

par(mfrow = c(1,1), mar = rep(0,4))
plot(tr, cex = 0.3)

write.tree(tr, "downloaded tree.phy")

##
## Run Models
##

# your_data$dam<-factor(your_data$dam, levels=levels(your_data$animal))
# invA <- inverseA(tr)
# ginverse=list(animal=invA, dam=invA)
summary(data$animal)

# Run full models

model2.2<- MCMCglmm( log_FMR ~ log_FMR_Mass + Phase + Latitude + log_Colony + Average_Brood,
                     random =~animal,
                     pedigree = tr,
                     pr = TRUE,
                     nodes="TIPS", scale=F,
                     data = data,  prior = prior, burnin = 300000/5, nitt = 1300000/5,
                     thin = 1000/5, verbose = T)

plot(model2.2)
summary(model2.2)
print(tr$tip.label)

model2.3<- MCMCglmm( log_FMR ~ log_FMR_Mass + Phase + Latitude + Average_Brood,
                     random = ~animal,
                     pedigree = tr,
                     pr = TRUE,
                     nodes="TIPS", scale=F,
                     data = data,  prior = prior,
                     burnin = 300000/5, nitt = 1300000/5,
                     thin = 1000/5, verbose = T)

plot(model2.3)
summary(model2.3)

model2.2$DIC
model2.3$DIC


# Run Group 3 models and incorporate random terms and ajust priors accordingly

model2.4<- MCMCglmm( log_FMR ~ log_FMR_Mass + Phase + Latitude + log_Colony + Average_Brood,
                     random = ~animal + Latin_Name,
                     pedigree = tr,
                     pr = TRUE,
                     nodes="TIPS", scale=F,
                     data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
                                                                G2=list(V=1,nu=0.001)),
                                                       R = list(V=1,nu=0.001)),
                     burnin = 300000/5, nitt = 1300000/5,
                     thin = 1000/5, verbose = T)

plot(model2.4)
summary(model2.4)

model2.5<- MCMCglmm( log_FMR ~ log_FMR_Mass + Phase + Latitude + Average_Brood,
                     random = ~animal + Latin_Name,
                     pedigree = tr,
                     pr = TRUE,
                     nodes="TIPS", scale=F,
                     data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
                                                                G2=list(V=1,nu=0.001)),
                                                       R = list(V=1,nu=0.001)),
                     burnin = 300000/5, nitt = 1300000/5,
                     thin = 1000/5, verbose = T)

plot(model2.5)
summary(model2.5)

model2.4$DIC
model2.5$DIC

# Run Group 4 models and incorporate random terms and ajust priors accordingly

model2.6 <- MCMCglmm( log_FMR ~ log_FMR_Mass + Phase + Latitude + log_Colony + Average_Brood,
                     random =~animal + Latin_Name + Site,
                     pedigree = tr,
                     pr = TRUE,
                     nodes="TIPS", scale=F,
                     data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
                                                                G2=list(V=1,nu=0.001),
                                                                G3=list(V=1,nu=0.001)),
                                                       R = list(V=1,nu=0.001)),
                     burnin = 300000/5, nitt = 1300000/5,
                     thin = 1000/5, verbose = T)

plot(model2.6)
summary(model2.6)

model2.7 <- MCMCglmm( log_FMR ~ log_FMR_Mass + Phase + Latitude + Average_Brood,
                      random =~animal + Latin_Name + Site,
                      pedigree = tr,
                      pr = TRUE,
                      nodes="TIPS", scale=F,
                      data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
                                                                 G2=list(V=1,nu=0.001),
                                                                 G3=list(V=1,nu=0.001)),
                                                        R = list(V=1,nu=0.001)),
                      burnin = 300000/5, nitt = 1300000/5,
                      thin = 1000/5, verbose = T)

plot(model2.7)
summary(model2.7)

model2.6$DIC
model2.7$DIC

#######################################################################################################

##
## Make Estimates
##

par(mar = c(4,4,2,0)+0.5)

#estimate for a brooding little auk
plot(model2.2$Sol[,"(Intercept)"] + model2.2$Sol[,"animal.alle alle"] + mean(data$Latitude)*model2.2$Sol[,"Latitude"]+
       mean(data$log_Colony)*model2.2$Sol[,"log_Colony"]+log10(150)*model2.2$Sol[,"log_FMR_Mass"] + model2.2$Sol[,"Phase2Brood"])
posterior.mode(model2.2$Sol[,"(Intercept)"] + model2.2$Sol[,"animal.alle alle"] + mean(data$Latitude)*model2.2$Sol[,"Latitude"]+
                 mean(data$log_Colony)*model2.2$Sol[,"log_Colony"]+log10(150)*model2.2$Sol[,"log_FMR_Mass"] + model2.2$Sol[,"Phase2Brood"])
HPDinterval(model2.2$Sol[,"(Intercept)"] + model2.2$Sol[,"animal.alle alle"] + mean(data$Latitude)*model2.2$Sol[,"Latitude"]+
              mean(data$log_Colony)*model2.2$Sol[,"log_Colony"]+log10(150)*model2.2$Sol[,"log_FMR_Mass"] + model2.2$Sol[,"Phase2Brood"])

#The red-legged kittiwake (Rissa brevirostris) 428g
posterior.mode(model2.2$Sol[,"(Intercept)"] + model2.2$Sol[,"animal.rissa brevirostris"] + mean(data$Latitude)*model2.2$Sol[,"Latitude"]+
                 mean(data$log_Colony)*model2.2$Sol[,"log_Colony"]+log10(428)*model2.2$Sol[,"log_FMR_Mass"] + model2.2$Sol[,"Phase2Brood"])
HPDinterval(model2.2$Sol[,"(Intercept)"] + model2.2$Sol[,"animal.rissa brevirostris"] + mean(data$Latitude)*model2.2$Sol[,"Latitude"]+
              mean(data$log_Colony)*model2.2$Sol[,"log_Colony"]+log10(428)*model2.2$Sol[,"log_FMR_Mass"] + model2.2$Sol[,"Phase2Brood"])

#Australian pelican (Pelecanus conspicillatus) 6120 g
posterior.mode(model2.2$Sol[,"(Intercept)"] + model2.2$Sol[,"animal.pelecanus conspicillatus"] + mean(data$Latitude)*model2.2$Sol[,"Latitude"]+
                 mean(data$log_Colony)*model2.2$Sol[,"log_Colony"]+log10(6120)*model2.2$Sol[,"log_FMR_Mass"] + model2.2$Sol[,"Phase2Brood"])
HPDinterval(model2.2$Sol[,"(Intercept)"] + model2.2$Sol[,"animal.pelecanus conspicillatus"] + mean(data$Latitude)*model2.2$Sol[,"Latitude"]+
              mean(data$log_Colony)*model2.2$Sol[,"log_Colony"]+log10(6120)*model2.2$Sol[,"log_FMR_Mass"] + model2.2$Sol[,"Phase2Brood"])

#Peruvian diving petrel (Pelecanoides garnotii) 207 g
posterior.mode(model2.2$Sol[,"(Intercept)"] + model2.2$Sol[,"animal.pelecanoides garnotii"] + mean(data$Latitude)*model2.2$Sol[,"Latitude"]+
                 mean(data$log_Colony)*model2.2$Sol[,"log_Colony"]+log10(207)*model2.2$Sol[,"log_FMR_Mass"] + model2.2$Sol[,"Phase2Brood"])
HPDinterval(model2.2$Sol[,"(Intercept)"] + model2.2$Sol[,"animal.pelecanoides garnotii"] + mean(data$Latitude)*model2.2$Sol[,"Latitude"]+
              mean(data$log_Colony)*model2.2$Sol[,"log_Colony"]+log10(207)*model2.2$Sol[,"log_FMR_Mass"] + model2.2$Sol[,"Phase2Brood"])

par(mfrow = c(1,1), mar = c(4,4,0,0)+0.5)
plot(FMR ~ FMR_Mass_g, data = data, type = "n", axes = F, log = "xy")

lines(x = range(subset(data, Phase == "1Incubation")$FMR_Mass_g),
      y = 10^(summary(model2.2)$solutions["(Intercept)","post.mean"]+
                mean(data$Latitude)*summary(model2.2)$solutions["Latitude","post.mean"]+
                mean(data$log_Colony)*summary(model2.2)$solutions["log_Colony","post.mean"]+
                log10(range(subset(data, Phase == "1Incubation")$FMR_Mass_g))*
                summary(model2.2)$solutions["log_FMR_Mass","post.mean"]),
      lwd = 1.5, col = pal.blue[5])
lines(x = range(subset(data, Phase == "2Brood")$FMR_Mass_g),
      y = 10^(summary(model2.2)$solutions["(Intercept)","post.mean"]+
                mean(data$Latitude)*summary(model2.2)$solutions["Latitude","post.mean"]+
                mean(data$log_Colony)*summary(model2.2)$solutions["log_Colony","post.mean"]+
                log10(range(subset(data, Phase == "2Brood")$FMR_Mass_g))*
                summary(model2.2)$solutions["log_FMR_Mass","post.mean"]+
                summary(model2.2)$solutions["Phase2Brood","post.mean"]),
      lwd = 1.5, col = pal.green[5], lty = 2)
lines(x = range(subset(data, Phase == "3Creche")$FMR_Mass_g),
      y = 10^(summary(model2.2)$solutions["(Intercept)","post.mean"]+
                mean(data$Latitude)*summary(model2.2)$solutions["Latitude","post.mean"]+
                mean(data$log_Colony)*summary(model2.2)$solutions["log_Colony","post.mean"]+
                log10(range(subset(data, Phase == "3Creche")$FMR_Mass_g))*
                summary(model2.2)$solutions["log_FMR_Mass","post.mean"]+
                summary(model2.2)$solutions["Phase3Creche","post.mean"]),
      lwd = 1.5, col = pal.red[5], lty = 3)


points(FMR ~ FMR_Mass_g, data = subset(data, Phase == "1Incubation"),
       pch = 21, bg = "white", col = pal.blue[5], cex = 1)
points(FMR ~ FMR_Mass_g, data = subset(data, Phase == "2Brood"),
       pch = 22, bg = "white", col = pal.green[5], cex = 1)
points(FMR ~ FMR_Mass_g, data = subset(data, Phase == "3Creche"),
       pch = 23, bg = "white", col = pal.red[5], cex = 1)

points(x = 428,
       y = 10^(posterior.mode(model2.2$Sol[,"(Intercept)"] + model2.2$Sol[,"animal.rissa brevirostris"] + mean(data$Latitude)*model2.2$Sol[,"Latitude"]+
                                mean(data$log_Colony)*model2.2$Sol[,"log_Colony"]+log10(428)*model2.2$Sol[,"log_FMR_Mass"] + model2.2$Sol[,"Phase2Brood"])),
       pch = 22, col = "white", bg = pal.green[5], cex = 1.5)
points(x = 6120,
       y = 10^(posterior.mode(model2.2$Sol[,"(Intercept)"] + model2.2$Sol[,"animal.pelecanus conspicillatus"] + mean(data$Latitude)*model2.2$Sol[,"Latitude"]+
                                mean(data$log_Colony)*model2.2$Sol[,"log_Colony"]+log10(6120)*model2.2$Sol[,"log_FMR_Mass"] + model2.2$Sol[,"Phase2Brood"])),
       pch = 22, col = "white", bg = pal.green[5], cex = 1.5)
points(x = 207,
       y = 10^(posterior.mode(model2.2$Sol[,"(Intercept)"] + model2.2$Sol[,"animal.pelecanoides garnotii"] + mean(data$Latitude)*model2.2$Sol[,"Latitude"]+
                                mean(data$log_Colony)*model2.2$Sol[,"log_Colony"]+log10(207)*model2.2$Sol[,"log_FMR_Mass"] + model2.2$Sol[,"Phase2Brood"])),
       pch = 22, col = "white", bg = pal.green[5], cex = 1.5)

axis(side = 1, las = 1, lwd = 1.5)
axis(side = 2, las = 1, lwd = 1.5)
box(lwd = 1.5)

#Extract phylogenetic heritability
model.for.variances <- model2.2

posterior.heritability.animal <- model.for.variances$VCV[, "animal"]/(model.for.variances$VCV[, "animal"] + model.for.variances$VCV[, "units"])
plot(posterior.heritability.animal)
HPDinterval(posterior.heritability.animal)
HPDinterval(model.for.variances$VCV[, "animal"])

