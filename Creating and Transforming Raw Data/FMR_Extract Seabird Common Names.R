
rm(list = ls(all = TRUE))

# install.packages("rredlist")

library(rredlist)
library(reshape2)

#

data <- read.csv("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/SeabirdSpecies.csv")

Sci_Name <- as.character(data$animal)
Sci_Name <- gsub("_", " ", Sci_Name, fixed=TRUE)

#

Names <- NULL

for (i in 1:nrow(data)){
  print(i)
  
  # i = 13

Name <- as.character(rl_common_names_(Sci_Name[i],
                        key = "2503fc54ee4a5fb0394b174d99901b4fd16bbc0485743d5fe9e36da28047072f",
                        parse = TRUE))

Names <- rbind(Names, data.frame(Name))
}

#

write.csv(species_choices, "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/Seabird_Sci_Common_Names.csv")

# 

# Once saved, open in excel, separate columns by "" marks and cut and paste common name only column into SeabirdSpecies.csv

#

# Load in (a) phylo tree and (b) seabird common names and families
# Merge 2 dataframes

tr <- read.tree("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/bird tree seabird tree.phy")
sptree <- makeNodeLabel(tr, method = "number", prefix = "node")
seabirds <- read.csv("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/SeabirdSpecies.csv")
species_choices <- as.data.frame(sort(sptree$tip.label))
names(species_choices)[1] <- "animal"
species_choices <- merge(species_choices, seabirds, by = "animal")
species_choices<-species_choices[,1:3]
rm(sptree, tr, seabirds)

#

species_choices$Sci_Name <- gsub("_", " ", species_choices$animal, fixed=TRUE)
species_choices$Common_Name <- paste("(",species_choices$Common_Name,")", sep = "")

species_choices$Both_Names <- paste(species_choices$Sci_Name, species_choices$Common_Name)
species_choices <- species_choices[,c(3,5,1)]
#

molten <- melt(species_choices, id.vars = "Family", measure.vars = "Both_Names")
hey <- dcast(molten, Family ~ variable, list)
