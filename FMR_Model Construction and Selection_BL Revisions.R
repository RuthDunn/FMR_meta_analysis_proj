
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
data <- read.xlsx("M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Data/Breeding FMR - Meta Analysis.xlsx", 1)

data$animal <- data$Match_name

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

################################################################################

# Select species to be included within phylogeny analysis

tr <- read.tree("M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Data/bird tree seabird tree.phy")

sptree <- makeNodeLabel(tr, method = "number", prefix = "node")
treeAinv<-inverseA(sptree, nodes="TIPS")$Ainv

################################################################################

##
## Run Models
##

# Incorporate random terms: phylogeny and species (Sci_Name) and colony (Short_Location)

# Adjust priors accordingly


# Full model:
# (Mass, Phase, Lat, colony size (log_colony_mass), Average Brood Size)

model1<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + log_colony_mass + Average_Brood,
                     random = ~animal + Sci_Name + Short_Location,
                     pedigree = sptree,
                     pr = TRUE,
                     nodes="TIPS", scale=F,
                     data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
                                                                G2=list(V=1,nu=0.001),
                                                                G3=list(V=1,nu=0.001)),
                                                       R = list(V=1,nu=0.001)),
                     burnin = 300000/5, nitt = 1300000/5,
                     thin = 1000/5, verbose = T)

save(model1, file = "M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Revised/model1.rda")

# load(file = "M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Revised/model2.rda")
summary(model1)
# plot(model2.1)


# Begin to remove non-significant terms:


# Edit:
# (Include: Mass, Phase, Lat, Pairs)
# (Remove: Average_Brood)

model2<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + log_colony_mass,
                     random = ~animal + Sci_Name + Short_Location,
                     pedigree = sptree,
                     pr = TRUE,
                     nodes="TIPS", scale=F,
                     data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
                                                                G2=list(V=1,nu=0.001),
                                                                G3=list(V=1,nu=0.001)),
                                                       R = list(V=1,nu=0.001)),
                     burnin = 300000/5, nitt = 1300000/5,
                     thin = 1000/5, verbose = T)
save(model2, file = "M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Revised/model2.rda")

# load(file = "M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Revised/model2.rda")
summary(model2)
# plot(model2)


# Seventh edit:
# (Include: Mass, Phase, Lat, Average Brood Size)
# (Remove: Pairs (log_colony_mass))

model3<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat + Average_Brood,
                   random = ~animal + Sci_Name + Short_Location,
                   pedigree = sptree,
                   pr = TRUE,
                   nodes="TIPS", scale=F,
                   data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
                                                              G2=list(V=1,nu=0.001),
                                                              G3=list(V=1,nu=0.001)),
                                                     R = list(V=1,nu=0.001)),
                   burnin = 300000/5, nitt = 1300000/5,
                   thin = 1000/5, verbose = T)
save(model3, file = "M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Revised/model3.rda")

# load(file = "M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Revised/model3.rda")
# summary(model3)
# plot(model3)


# Edit:
# (Include: Mass, Phase, Lat)
# (Remove: Pairs, Average_Brood)

model4<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat,
                     random = ~animal + Sci_Name + Short_Location,
                     pedigree = sptree,
                     pr = TRUE,
                     nodes="TIPS", scale=F,
                     data = data,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
                                                                G2=list(V=1,nu=0.001),
                                                                G3=list(V=1,nu=0.001)),
                                                       R = list(V=1,nu=0.001)),
                     burnin = 300000/5, nitt = 1300000/5,
                     thin = 1000/5, verbose = T)
save(model4, file = "M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Revised/model4.rda")

# load(file = "M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Revised/model4.rda")
# summary(model4)
# plot(model4)




################################################################################

# Compare model DICs:

model1$DIC
model2$DIC
model3$DIC
model4$DIC

################################################################################

summary(model1)
summary(model2)
summary(model3)
summary(model4)

################################################################################

phylo_h2.1 <- model1$VCV[,"animal"] / (model1$VCV[,"animal"] + model1$VCV[,"Sci_Name"] +
                                       model1$VCV[,"Short_Location"] + model1$VCV[,"units"])
summary(phylo_h2.1)

phylo_h2.2 <- model2$VCV[,"animal"] / (model2$VCV[,"animal"] + model2$VCV[,"Sci_Name"] +
                                       model2$VCV[,"Short_Location"] + model2$VCV[,"units"])
summary(phylo_h2.2)

phylo_h2.3 <- model3$VCV[,"animal"] / (model3$VCV[,"animal"] + model3$VCV[,"Sci_Name"] +
                                       model3$VCV[,"Short_Location"] + model3$VCV[,"units"])
summary(phylo_h2.3)

phylo_h2.4 <- model4$VCV[,"animal"] / (model4$VCV[,"animal"] + model4$VCV[,"Sci_Name"] +
                                       model4$VCV[,"Short_Location"] + model4$VCV[,"units"])
summary(phylo_h2.4)
