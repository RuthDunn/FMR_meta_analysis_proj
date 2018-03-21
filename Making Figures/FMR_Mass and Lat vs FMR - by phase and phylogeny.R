
rm(list = ls(all = TRUE))

# Colour blind friendly palette:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# install.packages(c("xlsx", "ape", "MCMCglmm"))


library(xlsx)
library(ape)
library(MCMCglmm)
library(ggplot2)
library(cowplot)

# prev.data <- read.csv("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/Breeding_FMR_Physiology.csv")
data <- read.xlsx("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/Breeding FMR - Meta Analysis.xlsx", 1)

data$animal <- data$Match_name
data$Colony_explicit <- as.numeric(data$Colony_explicit)
data$log_Colony <- log10(data$Colony_explicit)
data$Phase <- ordered(data$Phase, levels = c("Incubation", "Brood", "Creche"))

tr <- read.tree("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/bird tree seabird tree.phy")

sptree <- makeNodeLabel(tr, method = "number", prefix = "node")
treeAinv<-inverseA(sptree, nodes="TIPS")$Ainv

################################################################################

load(file = "C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/model2.5.rda")

# plot(model2.5)
# summary(model2.5)

################################################################################

##

plot(data$Lat, data$log_FMR, col = data$Phase)
abline(lm(data$log_FMR ~ data$Lat))

################################################################################

## 

mass.plot <- ggplot()

mass.plot <- mass.plot +
  geom_line(mapping = aes(x = range(subset(data, Phase == "Incubation")$Mass_g),
                          y = 10^(summary(model2.5)$solutions["(Intercept)","post.mean"]+
                                    mean(data$Lat)*summary(model2.5)$solutions["Lat","post.mean"]+
                                    log10(range(subset(data, Phase == "Incubation")$Mass_g))*
                                    summary(model2.5)$solutions["log_Mass","post.mean"]+
                                    summary(model2.5)$solutions["PhaseIncubation","post.mean"])), colour = "#E69F00") +
  
  geom_line(mapping = aes(x = range(subset(data, Phase == "Brood")$Mass_g),
                          y = 10^(summary(model2.5)$solutions["(Intercept)","post.mean"]+
                                    mean(data$Lat)*summary(model2.5)$solutions["Lat","post.mean"]+
                                    log10(range(subset(data, Phase == "Brood")$Mass_g))*
                                    summary(model2.5)$solutions["log_Mass","post.mean"])), colour = "#009E73") +
 
   geom_line(mapping = aes(x = range(subset(data, Phase == "Creche")$Mass_g),
                          y = 10^(summary(model2.5)$solutions["(Intercept)","post.mean"]+
                                    mean(data$Lat)*summary(model2.5)$solutions["Lat","post.mean"]+
                                    log10(range(subset(data, Phase == "Creche")$Mass_g))*
                                    summary(model2.5)$solutions["log_Mass","post.mean"]+
                                    summary(model2.5)$solutions["PhaseCreche","post.mean"])), colour = "#0072B2")

mass.plot <- mass.plot  +
  geom_point(mapping = aes(x = Mass_g, y = X.FMR, color = Phase, pch = Family), data = data) +
  scale_colour_manual(values=cbPalette[c(2, 4, 6)]) +
  scale_shape_manual(values = c(0:4, 6, 8, 15:17)) + 
  scale_x_log10("Mass (g)") +
  scale_y_log10("FMR (kJ/ day)")

hello <- theme_bw() %+replace% theme(panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank())

theme_set(hello)

mass.plot

# ggsave("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Plots/mass_fmr-phase_family.png",
#        width = 8, height = 5)

## 

lat.plot <- ggplot()

lat.plot <- lat.plot  +
  geom_point(mapping = aes(x = Lat, y = X.FMR, color = Phase, pch = Family), data = data) +
  scale_colour_manual(values=cbPalette[c(2, 4, 6)]) +
  scale_shape_manual(values = c(0:4, 6, 8, 15:17)) + 
  xlab("Latitude (North/ South)") +
  scale_y_log10(" ")

lat.plot <- lat.plot +
  geom_line(mapping = aes(x = range(subset(data, Phase == "Incubation")$Lat),
                          
                          y = 10^(summary(model2.5)$solutions["(Intercept)","post.mean"]+
                                    mean(data$log_Mass)*summary(model2.5)$solutions["log_Mass","post.mean"]+
                                    (range(subset(data, Phase == "Incubation")$Lat))*
                                    summary(model2.5)$solutions["Lat","post.mean"]+
                                    summary(model2.5)$solutions["PhaseIncubation","post.mean"])), colour = "#E69F00") +
  
  geom_line(mapping = aes(x = range(subset(data, Phase == "Brood")$Lat),
                          
                          y = 10^(summary(model2.5)$solutions["(Intercept)","post.mean"]+
                                    mean(data$log_Mass)*summary(model2.5)$solutions["log_Mass","post.mean"]+
                                    (range(subset(data, Phase == "Brood")$Lat))*
                                    summary(model2.5)$solutions["Lat","post.mean"])), colour = "#009E73") +
  
  geom_line(mapping = aes(x = range(subset(data, Phase == "Creche")$Lat),
                          
                          y = 10^(summary(model2.5)$solutions["(Intercept)","post.mean"]+
                                    mean(data$log_Mass)*summary(model2.5)$solutions["log_Mass","post.mean"]+
                                    (range(subset(data, Phase == "Creche")$Lat))*
                                    summary(model2.5)$solutions["Lat","post.mean"]+
                                    summary(model2.5)$solutions["PhaseCreche","post.mean"])), colour = "#0072B2")

lat.plot <- lat.plot + theme_set(hello)
# lat.plot + ylim(600, 950)

# lat.plot

##

prow <- plot_grid(mass.plot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0.4, 0), "cm")),
                  lat.plot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0.4, 0), "cm")),
                  labels = c("a)", "b)"))
legend <- get_legend(mass.plot)
p <- plot_grid(prow, legend, rel_widths = c(3, .7))
p

##

# ggsave("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Plots/mass_lat_fmr-phase_family.png",
#        width = 8, height = 5)

################################################################################
