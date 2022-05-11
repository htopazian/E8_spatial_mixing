# libraries
library(sf)
library(tidyverse)
library(patchwork)
library(terra)

# import admin 1 shapefiles from GADM
E8_SADC <- c('AGO','BWA','SWZ','MOZ','NAM','ZAF','ZMB','ZWE',
         'COM','COD','LSO','MDG','MWI','MUS','SYC','TZA')

SSA <- c('KEN','UGA','COG','SSD','ETH','SOM','CAF','CMR','GAB','GNQ','RWA','BDI','NGA')


# saves copy of gadm SpatVector as .rds in data folder
import_gadm <- function(ISO, level){
  geodata::gadm(country = ISO, level = level, path = './01_data/gadm', version="4.0")
}

map2(E8_SADC, 1, import_gadm)
map2(SSA, 0, import_gadm)

unpack_gadm <- function(file){
  object <- readRDS(file) # read in object
  object <- terra::vect(object) # unpack SpatVector
  st_as_sf(object) # transform to sf object
}

# create a list of all countries
files <- list.files(path = "./01_data/gadm", # Identify all .rds files in folder
                       pattern = "*.rds", full.names = TRUE)

# loop over each country
countries <- map_dfr(files, unpack_gadm)

st_crs(countries) # view CRS

E8 <- countries %>%
  filter(ID_0 %in% c('AGO','BWA','SWZ','MOZ','NAM','ZAF','ZMB','ZWE'))

SADC <- countries %>%
  filter(ID_0 %in% c('COM','COD','LSO','MDG','MWI','MUS','SYC','TZA'))

SSA <- countries %>%
  filter(ID_0 %in% SSA)


# find centroids in E8 admin1s
centroids <- sf::st_centroid(E8)


# plot
ggplot() +
  geom_sf(data = SSA, fill = "cornsilk2", color = "cornsilk3") +
  geom_sf(data = SADC, fill = "lightgreen", color = "cornsilk3") +
  geom_sf(data = E8, fill = "darkgreen") +
  geom_sf(data = centroids, color = 'red') +
  # geom_sf(data=E8, fill=NA, color="tan4", size=0.75) +
  theme_bw(base_size = 14) +
  scale_color_distiller(palette = 'Spectral') +
  scale_x_continuous(limits = c(10, 52)) +
  scale_y_continuous(limits = c(-35, 6)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "#daeff8", color = NA))

ggsave('./03_output/E8_map.pdf', width = 6, height = 6)
