# libraries
library(sf)
library(tidyverse)
library(patchwork)
library(terra)
library(geodata)        # gadm
library(wpgpDownloadR)  # worldpop devtools::install_github("wpgp/wpgpDownloadR")
library(malariaAtlas)   # malaria atlas project
library(raster)



# import GADM shapefiles -------------------------------------------------------

# define countries by ISO
# also downloadable at: https://geodata.ucdavis.edu/gadm/gadm4.0/pck/
E8_SADC <- c('AGO','BWA','SWZ','MOZ','NAM','ZAF','ZMB','ZWE',
         'COM','COD','LSO','MDG','MWI','MUS','SYC','TZA')

SSA <- c('KEN','UGA','COG','SSD','ETH','SOM','CAF','CMR','GAB','GNQ','RWA','BDI','NGA')

Senegambia <- c('SEN','GMB')

Senegambia_neighbors <- c('GNB','GIN','MRT','MLI')


# saves copy of gadm SpatVector as .rds in data folder
import_gadm <- function(ISO, level){
  geodata::gadm(country = ISO, level = level, path = './01_data/gadm', version="4.0")
}

map2(E8_SADC, 1, import_gadm)
map2(SSA, 0, import_gadm)
map2(Senegambia, 1, import_gadm)
map2(Senegambia_neighbors, 0, import_gadm)

# create a list of all countries
files <- list.files(path = "./01_data/gadm", # Identify all .rds files in folder
                       pattern = "*.rds", full.names = TRUE)

# unpack gadms
unpack_gadm <- function(file){
  object <- readRDS(file) # read in object
  object <- terra::vect(object) # unpack SpatVector
  st_as_sf(object) # transform to sf object
}

countries <- map_dfr(files, unpack_gadm) # loop over each country

saveRDS(countries, './03_output/countries.rds') # save
# countries <- readRDS('./03_output/countries.rds')

st_crs(countries) # view CRS

# partition data by country grouping
E8 <- countries %>%
  filter(ID_0 %in% c('AGO','BWA','SWZ','MOZ','NAM','ZAF','ZMB','ZWE'))

SADC <- countries %>%
  filter(ID_0 %in% c('COM','COD','LSO','MDG','MWI','MUS','SYC','TZA'))

SSA <- countries %>%
  filter(ID_0 %in% SSA)

Senegambia <- countries %>%
  filter(ID_0 %in% Senegambia)

Senegambia_neighbors <- countries %>%
  filter(ID_0 %in% Senegambia_neighbors)



# worldpop ---------------------------------------------------------------------
# https://github.com/wpgp/wpgpDownloadR

# download raster file and .csv of population counts
wpgpGetCountryDataset(ISO3 = 'SEN', covariate = 'ppp_2020', destDir = './01_data/worldpop')
wpgpGetCountryDataset(ISO3 = 'GMB', covariate = 'ppp_2020', destDir = './01_data/worldpop')

wpgpGetPOPTable(ISO3 = 'SEN', year = 2020, destDir = './01_data/worldpop')
wpgpGetPOPTable(ISO3 = 'GMB', year = 2020, destDir = './01_data/worldpop')

# import rasters and merge
SG_rast <- raster::merge(raster('./01_data/worldpop/sen_ppp_2020.tif'),
                         raster('./01_data/worldpop/gmb_ppp_2020.tif'))

# check res
res(SG_rast)

# downsize for ease of plotting
SG_rast <- aggregate(SG_rast, fact=3); res(SG_rast)

# extract population within each admin1
pop <- raster::extract(SG_rast, Senegambia, exact = T)

# sum for each polygon
pop <- unlist(lapply(pop, function(x) if (!is.null(x)) sum(x, na.rm = TRUE) else NA )) %>%
  as_tibble() %>% rename(wpop_pop = value)

saveRDS(pop, './03_output/Senegambia_population.rds') # save


# change to dataframe to use with ggplot()
SG_rast <- data.frame(rasterToPoints(SG_rast, spatial = T))



# MAP PfPR ---------------------------------------------------------------------
# https://malariaatlas.org/malaria-burden-data-download/
# Annual mean of PF Parasite Rate
pfpr <- read_csv('./01_data/map/00_PfPR_table_Global_admin1_2000-2019.csv') %>%
  filter(Year == 2019 & ISO %in% c('SEN', 'GMB')) %>%
  dplyr::select(ISO, Name_0, Name_1, Pop, PfPR_rmean) %>% rename(map_pop = Pop) %>%
  mutate(Name_1 = case_when(Name_1 == 'Thies' ~ 'Thiès',
                            Name_1 == 'Sedhiou' ~ 'Sédhiou',
                            Name_1 == 'Kedougou' ~ 'Kédougou',
                            Name_1 == 'West Coast' ~ 'Western',
                            Name_1 == 'Central River' ~ 'Maccarthy Island',
                            Name_1 == 'Kanifing Municipal Council' ~ 'Banjul',
                            TRUE ~ Name_1))

pfpr <- Senegambia %>% left_join(pfpr, by = c('ID_0' = 'ISO', 'NAME_1' = 'Name_1', 'COUNTRY' = 'Name_0'))

# calibrate pfprs to EIRs in odin
# https://github.com/mrc-ide/cali/blob/main/R/calibrate.R
calibrate <- function(target, target_tt, summary_function, tolerance, interval = c(0.01, 2000) / 365, ...){
  stats::uniroot(objective,
                 target = target,
                 target_tt = target_tt,
                 summary_function = summary_function,
                 tolerance = tolerance,
                 interval = interval,
                 ...)
}

objective <- function(x, target, target_tt, summary_function, tolerance){
  message("Trying EIR: ", signif(x, 3))

  raw_output <- ICDMM::run_model(model = "odin_model",
                                 init_EIR = x,
                                 init_ft = 0.4,
                                 time = target_tt)

  raw_output <- as_tibble(raw_output)

  target_variable <- raw_output$prev2to10
  difference <- (target_variable[target_tt] - target)
  message("Current difference: ", paste(signif(difference, 3), collapse = " "))
  # Adjust for specified tolerance
  difference[abs(difference) < tolerance] <- 0
  difference <- sum(difference)
  return(difference)
}

summary_pfpr_2_10 <- function(x){
  prev_2_10 <- x$prev2to10
  return(prev_2_10)
}


PR_EIR <- function(target, target_tt){
  set.seed(123)
  out <- calibrate(target = target,
                   target_tt = target_tt,
                   summary_function = summary_pfpr_2_10,
                   tolerance = 0.02,
                   interval = c(.0001, 300))

  root <- as_tibble(out) %>% dplyr::select(root)

  root

}

# run over prevalence values for each admin1
EIR <- map2_dfr(pfpr$PfPR_rmean, rep(365, nrow(pfpr)), PR_EIR)


# aggregate pfpr and EIR estimates
pfpr_eir <- cbind(pfpr, EIR)

saveRDS(pfpr_eir, './03_output/Senegambia_pfpr_eir.rds') # save


# save master dataset of variables
master <- cbind(pfpr_eir, pop) %>%
  rename(EIR = root) %>%
  dplyr::select(ID_0, COUNTRY, NAME_1, VARNAME_1, map_pop, wpop_pop,
                PfPR_rmean, EIR, geometry)

saveRDS(master, './03_output/Senegambia_master.rds') # save


# E8 map -----------------------------------------------------------------------
# find centroids in E8 admin1s
centroids <- sf::st_centroid(E8)

# plot
ggplot() +
  geom_sf(data = SSA, fill = "cornsilk2", color = "cornsilk3") +
  geom_sf(data = SADC, fill = "lightgreen", color = "cornsilk3") +
  geom_sf(data = E8, fill = "darkgreen") +
  geom_sf(data = centroids, color = 'red', shape = 4, size = 2) +
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



# Senegambia map ---------------------------------------------------------------
# find centroids in E8 admin1s
centroids <- sf::st_centroid(Senegambia)
saveRDS(centroids, './03_output/Senegambia_centroids.rds')

# plot
ggplot() +
  geom_sf(data = Senegambia_neighbors, fill = "cornsilk2", color = "cornsilk3") +
  geom_sf(data = Senegambia[Senegambia$ID_0 == 'GMB',], fill = "#72B000", color = "cornsilk3") +
  geom_sf(data = Senegambia[Senegambia$ID_0 == 'SEN',], fill = "#72B000", color = "cornsilk3") +
  geom_raster(data = SG_rast, aes(x = x, y = y, alpha = layer), fill = "darkgreen", show.legend = F) +
  geom_sf(data = centroids, color = 'red', shape = 4, size = 2) +
  theme_bw(base_size = 14) +
  scale_alpha_continuous(range = c(0, 300), breaks = c(0, 0.1, 0.2, 0.5, 1, 300)) +
  scale_x_continuous(limits = c(-17.8, -11)) +
  scale_y_continuous(limits = c(12, 17)) +
  labs(x = '', y = '') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "#daeff8", color = NA))

ggsave('./03_output/Senegambia.pdf', width = 6, height = 6)


