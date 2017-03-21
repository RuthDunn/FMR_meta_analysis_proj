rm(list = ls(all = TRUE))

data <- read.csv("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/Breeding_FMR_Physiology_2.csv")
# data <- read.csv("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/")


data$animal <- as.character(data$animal)
data$Latitude <- abs(data$Latitude)
data$Colony_explicit <- as.numeric(data$Colony_explicit)
data$log_Colony <- log10(data$Colony_explicit)

library(rotl)
# library(nlme)
# library(lme4)
# library(MCMCglmm)
# library(ape)
library(car)
# library(RODBC)
# library(rncl)
library(phytools)
library(ggplot2)

# install.packages("SparseM")

# qplot(log_FMR_Mass, log_FMR, data = data,
#       ylab = "Log Field Metabolic Rate", xlab = "Log Mass (kg)",
#       color = Latitude,
#       facets = colvar~.Phase)
# 
# # Linear model of FMR and variables
# m2.1 <- lm(log_FMR ~ log_FMR_Mass + Phase +  Latitude + log_Colony + Average_Brood, 
#            data = data, na.action = na.omit)
# 
# # Plots of FMR against each variables, taking into account the effect of the other variables
# avPlots(m2.1)

###
### Select species to be included within phylogeny analysis
###

# birdtree <- read.nexus("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/BirdTreePhylogeny/2102.tre")
# plot(birdtree)
# cophenetic(tree)

# Red-legged kittiwake (Rissa brevirostris)
# Australian pelican (Pelecanus conspicillatus)
# Peruvian diving petrel (Pelecanoides garnotii)
# extra.species <- c("Rissa brevirostris", "Pelecanus conspicillatus", "Pelecanoides garnotii")
# resolved_names <- tnrs_match_names(c(data$animal,extra.species))

# Set this to TRUE to download the tree from ROTL
download.tree <- TRUE

# Use seabirds ott_ids (Open Tree Taxonomy identifier) to match and extract seabirds from the Open Tree tree
# seabirds <- read.csv("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/SeabirdSpecies.csv")
# seabirds$animal <- as.character(seabirds$animal)

# Try to load all seabird species
# resolved_names <- tnrs_match_names(data$animal)
# resolved_names <- tnrs_match_names(seabirds$animal)
# Error: "Queries containing more than 250 strings are not supported.
# You may submit multiplesmaller queries to avoid this limit"

# Try mutliple smaller queries:
# resolved_names1 <- tnrs_match_names(head(seabirds$animal,250))
# resolved_names2 <- tnrs_match_names(tail(seabirds$animal,101))
# resolved_names3 <- tnrs_match_names(data$animal)
# resolved_names <- rbind(resolved_names1, resolved_names2)
# rm(resolved_names1, resolved_names2,resolved_names3)
# 
# # Get the tree with just those tips:
# tr <- tol_induced_subtree(ott_ids=ott_id(resolved_names))
# 
# 
# # Open Tree trees have no branch lengths but MCMCglmm accepts topology as input
# taxon_map <- structure (resolved_names$search_string, names = resolved_names$unique_name)
# 
# 
# 
# # Tree contains ID numbers which don't therefore match our data
# # Use function to remove extra info at label tips
# # Then use taxon map to replace the tip labels with sp. names from the data
# otl_tips <- strip_ott_ids(tr$tip.label, remove_underscores=TRUE)
# tr$tip.label <- taxon_map[ otl_tips ]
# # Remove node labels by setting node.label to null
# tr$node.label <- NULL


###########

# Included sp. only
resolved_names3 <- tnrs_match_names(data$animal)

tr_sponly <- tol_induced_subtree(ott_ids=ott_id(resolved_names3))

taxon_map <- structure (resolved_names3$search_string, names = resolved_names3$unique_name)
otl_tips <- strip_ott_ids(tr_sponly$tip.label, remove_underscores=TRUE)
tr_sponly$tip.label <- taxon_map[ otl_tips ]
# Remove node labels by setting node.label to null
tr_sponly$node.label <- NULL



# write.tree(tr, "downloaded tree.phy")

# par(mfrow = c(1,1), mar = rep(0,4))
# plot(tr, cex = 0.3, type = "fan")

#
## Plotting
#

# Colours: http://paletton.com/#uid=63A0u0kl1Wx1x+IcEXDsUWkWEVB
# Primary color:
# pal.blue  <- ("#5FA8F3")
# Secondary color:
pal.green   <- ("#EBEBFF")
# Complement color:
pal.orange    <- ("#FC9428")

library(ggplot2)
library(phyloseq)
library(ggtree)
library(colorspace)
library(Cairo)
library(png)
library(RCurl)
library(EBImage)
library(Matrix)
library(rncl)
library(nlme)
library(cowplot)

# p <- ggtree(tr2, layout = "circular", branch.length = "none") +
#   geom_tippoint(size = 3, alpha = 1/4, aes(color = group)) +
#   geom_tiplab2(aes(subset = group, angle = angle), size = 1.5) +
#   scale_colour_manual(values = c(pal.green, rainbow(24)))
# 
# p <- p +
#   theme(
#     rect = element_rect(fill = "transparent") # bg of the panel
#   )
# p

# ggsave("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Phylo_Plot.png", bg = "transparent")




tr_sponly <- groupClade(tr_sponly, node = c(47, 54, 56, 59, 64, 63, 67, 69, 74, 75))

p <- ggtree(tr_sponly, color = pal.green, aes(size = 4), layout = "circular", branch.length = "none") +
  geom_tippoint(size = 6, alpha = 1/4, aes(color = group)) +
#  geom_tiplab2(aes(angle = angle), size = 3.8) +
#  geom_text2(aes(subset=!isTip, label=node)) +
  ggplot2::xlim(0, 10)


p <- p +
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent") # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
  #  , panel_border(remove = TRUE)
    , legend.background = element_rect(fill = "transparent") # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

p


# ggsave("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Phylo_Plot_names.png", bg = "transparent",
#     width = 30, height = 30, units= "cm")




