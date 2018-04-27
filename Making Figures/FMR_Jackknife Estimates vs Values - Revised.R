
rm(list = ls(all = TRUE))

################################################################################

library(ggplot2)
library(gplots)

All.estimates <- read.csv("M:/Documents/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Models/Jackknifing/Jackknife_Estimates-Post_Review.csv")

All.estimates$Phase <- ordered(All.estimates$Phase, levels = c("Incubation", "Brood", "Creche"))

# Simple plots:

# plot(All.estimates$estimate, Incubation$log_FMR)
# 
# plotCI(x = All.estimates$log_FMR, y = All.estimates$estimate,
#        ui = All.estimates$upper,
#        li = All.estimates$lower,
#        xlab = "log FMR Extracted Value",
#        ylab = "log FMR Estimates", gap = FALSE)
# abline(0,1, lty = 2)

# ggplot:

# Colour blind friendly palette:

hello <- theme_bw() %+replace% theme(panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank())

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot <- ggplot(All.estimates, aes(x = log_FMR, y = estimate))  +
  geom_abline(linetype = 2) +
  geom_point(aes(x = log_FMR, y = estimate, color = Phase)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.01, color = Phase)) +
  scale_colour_manual(values=cbPalette[c(7, 4, 6)]) +
  # scale_shape_manual(values = c(0:4, 6, 8, 15:17)) + 
  xlab("True log FMR Values (kJ/ day)")+
  ylab("Jackknife Estimate log FMR Values (kJ/ day)")

theme_set(hello)

plot

ggsave("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Plots/jackknife_estimates2.png",
       width = 8, height = 5)
