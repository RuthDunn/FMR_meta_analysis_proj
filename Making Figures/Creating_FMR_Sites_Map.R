# install.packages("ggmap")

library(rgdal)
library(ggmap)
library(maptools)

data <- read.csv("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Data/Breeding_FMR_Physiology_2.csv")
world <- readShapeSpatial("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Mapping/TM_WORLD_BORDERS-0.3.shp",
                          proj4string = CRS("+proj=longlat +datum=WGS84"))

world.data <- fortify(world)

get_map(location = "world")
ggmap(location = "world", data = data)


mp <- NULL
mapWorld <- borders("world", colour="transparent", fill="#EBEBFF") # create a layer of borders
mp <- ggplot() +   mapWorld

#Now Layer the cities on top
mp <- mp + geom_point(aes(x = Longitude, y = Latitude)
                      , color = Family, size=3,
                      data = data) 


mp <- mp +
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent") # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    #  , panel_border(remove = TRUE)
    , legend.background = element_rect(fill = "transparent") # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

mp


# ggsave("C:/Users/RuthDunn/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/Outputs/Sites_Map.png", bg = "transparent",
#     width = 30, height = 30, units= "cm")


