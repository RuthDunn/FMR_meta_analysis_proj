
rm(list = ls(all = TRUE))

# Sort out packages:

# install.packages(c("rJava", "xlsx", "ape", "MCMCglmm", "car", "openxlsx"))

library(rJava)
library(xlsx)
library(ape)
library(MCMCglmm)
library(car)
library(openxlsx)

# Load data:

data <- read.xlsx("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/Breeding FMR - Meta Analysis.xlsx", 1)

data$animal <- data$Match_name
data$Colony_explicit <- as.numeric(data$Colony_explicit)
data$log_Colony <- log10(data$Colony_explicit)

tr <- read.tree("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/bird tree seabird tree.phy")

sptree <- makeNodeLabel(tr, method = "number", prefix = "node")
treeAinv<-inverseA(sptree, nodes="TIPS")$Ainv

#

# Loop Model
# Include: Mass, Phase, Lat)

for (i in 21:nrow(data)){
  print(i)
  
  # i = 15

  data.drop = data[-i,]
  
mod<- MCMCglmm( log_FMR ~ log_Mass + Phase + Lat,
                     random = ~animal + Sci_Name,
                     pedigree = sptree,
                     pr = TRUE,
                     nodes="TIPS", scale=F,
                     data = data.drop,  prior = prior<-list(G = list(G1=list(V=1,nu=0.001),
                                                                G2=list(V=1,nu=0.001)),
                                                       R = list(V=1,nu=0.001)),
                     burnin = 300000/5, nitt = 1300000/5,
                     thin = 1000/500, verbose = T)

# save(mod, file = paste("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Jackknifing/",
#                        "mod_", i, ".rda", sep = ""))

save(mod, file = paste("M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Jackknifing/",
                       "mod_", i, ".rda", sep = ""))

# save(mod, file = paste("F:/Jackknifing/",
#                        "mod_", i, ".rda", sep = ""))

}
