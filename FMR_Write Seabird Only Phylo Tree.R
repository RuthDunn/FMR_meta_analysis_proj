
rm(list = ls(all = TRUE))


## Write seabird phylo tree

library(ape)

# Load 346 seabird species:
# (minus ducks, grebes, loons and phalaropes = 319 sp)

seabirds <- read.csv("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/SeabirdSpecies.csv")
seabirds$animal <- as.character(seabirds$animal)

# Load Ericson backbone tree from Jtz, Thomas, Joy, Harmann & Moores 2012
# (downloaded from birdtree.org)

MyTree <- read.tree((file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/BirdTree/AllBirdsEricson1.tre"))
sptree <- makeNodeLabel(MyTree[[1]], method = "number", prefix = "node")

# Select seabirds only

seabird.tree<-drop.tip(sptree, setdiff(sptree$tip.label, seabirds$animal))

# Write tree

write.tree(seabird.tree, "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/bird tree seabird tree.phy")