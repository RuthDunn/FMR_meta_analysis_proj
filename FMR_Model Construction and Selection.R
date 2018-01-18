
rm(list = ls(all = TRUE))

# Sort out packages:

# install.packages(c("xlsx", "ape", "MCMCglmm"))

library(rJava)
library(xlsx)
library(ape)
library(MCMCglmm)
library(car)

# library(nlme)
# library(lme4)
# library(RODBC)
# library(rncl)
# library(tidyr)
# library(phytools)
# library(MASS)

# Load data:

# prev.data <- read.csv("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/Breeding_FMR_Physiology.csv")
data <- read.xlsx("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/Breeding FMR - Meta Analysis.xlsx", 1)

data$animal <- data$Match_name
data$Colony_explicit <- as.numeric(data$Colony_explicit)
data$log_Colony <- log10(data$Colony_explicit)

################################################################################

# Basic model set up:

      # model_simple <- MCMCglmm(phen ~ cofactor,
          # random = ~phylo,
          # family = "gaussian", ginverse = list(phylo = inv.phylo$Ainv), prior=prior,
          # data = data, nitt = 5000000, burnin = 1000, thin = 500)

# These are "flat", "noninformative" priors:
prior <- list(G = list(G1 = list(V = 1, nu = 0.001)),    # Random effect
            R = list(V = 1, nu = 0.001))                 # Fixed effect

# Nitt = Total number of runs > 10,000
# BurnIn (running time before) = Recommended to be 10% of iterations and > 1000
# Thin = Sampling frequency

# Model without phylogeny effect:

# model0.1<- MCMCglmm(log_FMR ~ log_Mass + Phase +  Latitude + log_Colony + Average_Brood,
#                      random =~animal,
#                      data = data,  prior = prior,
#                      burnin = 10000, nitt = 100000, thin = 200, verbose = T)

# summary(model0.1)
# plot(model0.1)

# ^ Left hand plots should show no autocorrelation
# ^ Right hand plots should show normal distributions

# rm(model0.1)


# Linear model of FMR and variables:

# m0.1 <- lm(log_FMR ~ log_Mass + Phase +  Latitude + log_Colony + Average_Brood,
#            data = data, na.action = na.omit)

# Plots of FMR against each variables, taking into account the effect of the other variables:

# avPlots(m0.1)
# rm(m0.1)

################################################################################

# Select species to be included within phylogeny analysis

tr <- read.tree("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/bird tree seabird tree.phy")

sptree <- makeNodeLabel(tr, method = "number", prefix = "node")
treeAinv<-inverseA(sptree, nodes="TIPS")$Ainv

################################################################################

##
## Run Models
##

# Random term: phylogeny only

# Full model:
# (Mass, Phase, Lat, Long, Pairs (log_Colony), Average Brood Size)

# model1.1<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + Long + log_Colony + Average_Brood,
#                      random =~animal,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior, burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model1.1, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.1.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.1.rda")
summary(model1.1)
# plot(model1.1)


# Begin to remove non-significant terms:
# First edit:
# (Include: Mass, Phase, Lat, Long, Average Brood Size)
# (Remove: Pairs (log_Colony))

# model1.2<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + Long + Average_Brood,
#                      random = ~animal,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior,
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model1.2, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.2.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.2.rda")
summary(model1.2)
# plot(model1.2)


# Second edit:
# (Include: Mass, Phase, Lat, Long)
# (Remove: Pairs and Average_Brood)

# model1.3<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + Long,
#                      random = ~animal,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior,
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model1.3, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.3.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.3.rda")
summary(model1.3)
# plot(model1.3)


# Third edit:
# (Include: Mass, Phase, Lat, Average Brood)
# (Remove: Pairs, Long)

# model1.4<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + Average_Brood,
#                      random = ~animal,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior,
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model1.4, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.4.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.4.rda")
summary(model1.4)
# plot(model1.4)


# Fourth edit:
# (Include: Mass, Phase, Lat)
# (Remove: Pairs, Long, Average_Brood)

# model1.5<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat,
#                      random = ~animal,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior,
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model1.5, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.5.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.5.rda")
summary(model1.5) 
plot(model1.5)


# Fifth edit:
# (Include: Mass, Phase, Lat, Pairs)
# (Remove: Long, Average Brood) 

# model1.6<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + log_Colony,
#                      random = ~animal,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior,
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model1.6, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.6.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.6.rda")
summary(model1.6)
# plot(model1.6)


# Sixth edit:
# (Include: Mass, Phase, Lat, Pairs, Average Brood)
# (Remove: Long)

# model1.7<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + log_Colony + Average_Brood,
#                      random = ~animal,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior,
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model1.7, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.7.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.7.rda")
summary(model1.7)
# plot(model1.7)


# Seventh edit:
# (Include: Mass, Phase, Lat, Pairs, Long)
# (Remove: Average_Brood)

# model1.8<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + log_Colony + Long,
#                      random = ~animal,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior,
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model1.8, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.8.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model1.8.rda")
summary(model1.8)
# plot(model1.8)


# Compare models using DIC

model1.1$DIC
model1.2$DIC
model1.3$DIC
model1.4$DIC
model1.5$DIC
model1.6$DIC
model1.7$DIC
model1.8$DIC

# rm(model1.1, model1.2, model1.3, model1.4, model1.5, model1.6, model1.7, model1.8)

################################################################################

# Run "Group 3" models
# Incorporate random terms: phylogeny and species (Sci_Name)
# Adjust priors accordingly


# Full model:
# (Mass, Phase, Lat, Long, Pairs (log_Colony), Average Brood Size)

# model2.1<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + log_Colony + Average_Brood + Long,
#                      random = ~animal + Sci_Name,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
#                                                                 G2=list(V=1,nu=0.001)),
#                                                        R = list(V=1,nu=0.001)),
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model2.1, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.1.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.1.rda")
summary(model2.1)
# plot(model2.1)


# Begin to remove non-significant terms:
# First edit:
# (Include: Mass, Phase, Lat, Pairs, Long)
# (Remove: Average_Brood)

# model2.2<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + Long + log_Colony,
#                      random = ~animal + Sci_Name,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
#                                                                 G2=list(V=1,nu=0.001)),
#                                                        R = list(V=1,nu=0.001)),
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model2.2, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.2.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.2.rda")
summary(model2.2)
# plot(model2.2)


# Second edit:
# (Include: Mass, Phase, Lat, Pairs)
# (Remove: Long, Average Brood)

# model2.3<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + log_Colony,
#                      random = ~animal + Sci_Name,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
#                                                                 G2=list(V=1,nu=0.001)),
#                                                        R = list(V=1,nu=0.001)),
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model2.3, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.3.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.3.rda")
summary(model2.3)
# plot(model2.3)


# Third edit:
# (Include: Mass, Phase, Lat, Long)
# (Remove: Pairs (log_Colony) and Average_Brood)

# model2.4<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + Long,
#                      random = ~animal + Sci_Name,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
#                                                                 G2=list(V=1,nu=0.001)),
#                                                        R = list(V=1,nu=0.001)),
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model2.4, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.4.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.4.rda")
summary(model2.4)
# plot(model2.4)


# Fourth edit:
# (Include: Mass, Phase, Lat)
# (Remove: Pairs, Long, Average_Brood)

# model2.5<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat,
#                      random = ~animal + Sci_Name,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
#                                                                 G2=list(V=1,nu=0.001)),
#                                                        R = list(V=1,nu=0.001)),
#                      burnin = 300000/50, nitt = 1300000/50,
#                      thin = 1000/50, verbose = T)
# save(model2.5, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.5.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.5.rda")
summary(model2.5)
# plot(model2.5)


# Fifth edit:
# (Include: Mass, Phase, Lat, Average Brood)
# (Remove: Pairs, Long)

# model2.6<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + Average_Brood,
#                      random = ~animal + Sci_Name,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
#                                                                 G2=list(V=1,nu=0.001)),
#                                                        R = list(V=1,nu=0.001)),
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model2.6, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.6.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.6.rda")
summary(model2.6)
# plot(model2.6)


# Sixth edit:
# (Include: Mass, Phase, Lat, Pairs, Average Brood)
# (Remove: Long)

# model2.7<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + Average_Brood + log_Colony,
#                      random = ~animal + Sci_Name,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
#                                                                 G2=list(V=1,nu=0.001)),
#                                                        R = list(V=1,nu=0.001)),
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model2.7, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.7.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.7.rda")
summary(model2.7) 
# plot(model2.7)


# Seventh edit:
# (Include: Mass, Phase, Lat, Long, Average Brood Size)
# (Remove: Pairs (log_Colony))

# model2.8<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + Average_Brood + Long,
#                      random = ~animal + Sci_Name,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
#                                                                 G2=list(V=1,nu=0.001)),
#                                                        R = list(V=1,nu=0.001)),
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model2.8, file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.8.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.8.rda")
summary(model2.8)
# plot(model2.8)

# Eighth edit:
# (Include: Mass, Lat)
# (Remove: Phase, Pairs, Long, Average_Brood)

# model2.9<- MCMCglmm( log_FMR ~ log_Mass + Lat,
#                      random = ~animal + Sci_Name,
#                      pedigree = sptree,
#                      pr = TRUE,
#                      nodes="TIPS", scale=F,
#                      data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
#                                                                 G2=list(V=1,nu=0.001)),
#                                                        R = list(V=1,nu=0.001)),
#                      burnin = 300000/5, nitt = 1300000/5,
#                      thin = 1000/5, verbose = T)
# save(model2.9, file = "C:/Users/Ruth/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.9.rda")

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.9.rda")
summary(model2.9)
# plot(model2.9)

# Compare model DICs:

model2.1$DIC
model2.2$DIC
model2.3$DIC     # ****
model2.4$DIC
model2.5$DIC
model2.6$DIC
model2.7$DIC
model2.8$DIC
model2.9$DIC

################################################################################
