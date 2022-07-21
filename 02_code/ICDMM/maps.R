# creating a master file of parameters for Senegambia admin1s

# libraries
library(sf)
library(tidyverse)
library(patchwork)
library(terra)
library(geodata)        # gadm
library(wpgpDownloadR)  # worldpop devtools::install_github("wpgp/wpgpDownloadR")
library(malariaAtlas)   # malaria atlas project
library(raster)
library(umbrella)       # devtools::install_github('https://github.com/mrc-ide/umbrella')


# GADM shapefiles --------------------------------------------------------------

# define countries by ISO
# also downloadable at: https://geodata.ucdavis.edu/gadm/gadm4.0/pck/
E8_SADC <- c('AGO','BWA','SWZ','MOZ','NAM','ZAF','ZMB','ZWE',
         'COM','COD','LSO','MDG','MWI','MUS','SYC','TZA')

Senegambia <- c('SEN','GMB')

Senegambia_neighbors <- c('GNB','GIN','MRT','MLI')

# other countries in SSA
SSA <- c('KEN','UGA','COG','SSD','ETH','SOM','CAF','CMR','GAB','GNQ','RWA','BDI','NGA')


# saves copy of gadm SpatVector as .rds in data folder
import_gadm <- function(ISO, level){
  geodata::gadm(country = ISO, level = level, path = './01_data/gadm', version="4.0")
}

map2(E8_SADC, 1, import_gadm)              # admin1
map2(SSA, 0, import_gadm)                  # admin0
map2(Senegambia, 1, import_gadm)           # admin1
map2(Senegambia_neighbors, 0, import_gadm) # admin0

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

# downsize resolution for ease of plotting
SG_rast <- aggregate(SG_rast, fact=3); res(SG_rast)

# extract population within each admin1
pop <- raster::extract(SG_rast, Senegambia, exact = T)

# sum for each polygon
pop_2020 <- unlist(lapply(pop, function(x) if (!is.null(x)) sum(x, na.rm = TRUE) else NA )) %>%
  as_tibble() %>%
  rename(wpop_pop = value)


# weighted polygon centroids
# convert polygons to raster layer
z <- rasterize(Senegambia, SG_rast)
# compute weighted x and y coordinates within each rasterized region
xx <- zonal(init(SG_rast, v = "x") * SG_rast, z) / zonal(SG_rast, z) # default = mean
yy <- zonal(init(SG_rast, v = "y") * SG_rast, z) / zonal(SG_rast, z) # default = mean

# combine results in a matrix
weighted_centroids <- tibble(xx[, 2], yy[, 2])
head(weighted_centroids)


# change to dataframe to use with ggplot()
SG_rast <- data.frame(rasterToPoints(SG_rast, spatial = T))



# MAP ITN ----------------------------------------------------------------------
ITN_2020 <- read_csv('./01_data/map/00_ITN_table_Africa_admin1_2000-2020.csv') %>%
  filter(ISO %in% c('SEN', 'GMB') & Year == 2020) %>%
  dplyr::select(ISO, Name_1, Mean_ITN_coverage_rate) %>%
  mutate(Name_1 = case_when(Name_1 == 'Thies' ~ 'Thiès',
                            Name_1 == 'Sedhiou' ~ 'Sédhiou',
                            Name_1 == 'Kedougou' ~ 'Kédougou',
                            Name_1 == 'West Coast' ~ 'Western',
                            Name_1 == 'Central River' ~ 'Maccarthy Island',
                            Name_1 == 'Kanifing Municipal Council' ~ 'Banjul',
                            TRUE ~ Name_1))


# MAP PfPR ---------------------------------------------------------------------
# https://malariaatlas.org/malaria-burden-data-download/
# Annual mean of PF Parasite Rate
pfpr <- read_csv('./01_data/map/00_PfPR_table_Global_admin1_2000-2019.csv') %>%
  filter(Year == 2020 & ISO %in% c('SEN', 'GMB')) %>%
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
# not incorporating ITNs yet ...
source('./02_code/ICDMM/calibrate.R')

# run over prevalence values for each admin1
EIR <- map2_dfr(pfpr$PfPR_rmean, rep(365, nrow(pfpr)), PR_EIR)


# aggregate pfpr and EIR estimates
pfpr_eir <- cbind(pfpr, EIR)

saveRDS(pfpr_eir, './03_output/Senegambia_pfpr_eir.rds') # save
# pfpr_eir <- readRDS('./03_output/Senegambia_pfpr_eir.rds')


# Centroids --------------------------------------------------------------------

# find centroids in E8 admin1s
centroids <- sf::st_centroid(Senegambia)
saveRDS(centroids, './03_output/Senegambia_centroids.rds')

weighted_centroids <- st_as_sf(weighted_centroids, coords = c('xx[, 2]', 'yy[, 2]'))
st_crs(weighted_centroids) <- st_crs(centroids)

saveRDS(weighted_centroids, './03_output/Senegambia_weighted_centroids.rds')

centroids_combo <- rbind(
  centroids %>% dplyr::select(geometry) %>% mutate(type = 'geographic centroids'),
  weighted_centroids %>% mutate(type = 'weighted centroids')
)



# Seasonality ------------------------------------------------------------------

# specify time-range
start_date <- "2020-01-01"
end_date <- "2020-12-31"

# extract
daily_rain_raw <- umbrella::pull_daily_rainfall(sf = Senegambia,
                                      start_date = start_date,
                                      end_date = end_date)

# process the raw data
daily_rain  <- daily_rain_raw %>%
  # Convert to long and format
  pivot_longer(-c(ID_0:HASC_2),
               names_to = "date",
               values_to = "rainfall",
               names_prefix = "X",
               names_pattern = "(.*)_precipitation") %>%
  mutate(date = as.Date(as.character(readr::parse_number(.data$date)), format = "%Y%m%d"),
         year = lubridate::year(.data$date),
         day_of_year = lubridate::yday(.data$date),
         t = lubridate::yday(.data$date) / 365,
         rainfall = as.numeric(rainfall)) %>%
  # Replace missing data with 0
  replace_na(replace = list(rainfall = 0)) %>%
  # Remove any leap year addtional days
  dplyr::filter(.data$day_of_year < 366)

rain <- daily_rain %>%
  group_by(COUNTRY, NAME_1) %>%
  summarise(
    data = list(data.frame(date, day_of_year, t, rainfall)),
    model = list(umbrella::fit_fourier(data[[1]]$rainfall, data[[1]]$t)),
    profile = list(umbrella::fourier_predict(model[[1]]$coefficients, t = 1:365 / 365, floor = model[[1]]$floor))
  )

# extract model estimates
seasonality <- map_dfr(.x = 1:nrow(rain),
                       function(x){
                          tibble(list(rain[["model"]][[x]][["estimate"]]))
                         })

seasonality <- cbind(Senegambia$COUNTRY, Senegambia$NAME_1, seasonality)

colnames(seasonality) <- c('COUNTRY', 'NAME_1', 'seasonality')


# Save master dataset ----------------------------------------------------------
master <- cbind(pfpr_eir, pop_2020) %>%
  rename(EIR = root) %>%
  dplyr::select(ID_0, COUNTRY, NAME_1, VARNAME_1, map_pop, wpop_pop,
                PfPR_rmean, EIR, geometry)
master <- cbind(master,
                weighted_centroids %>% rename(w_cent_geometry = geometry),
                centroids %>% dplyr::select(geometry) %>% rename(cent_geometry = geometry))

master <- master %>% left_join(ITN_2020, by = c('ID_0' = 'ISO', 'NAME_1' = 'Name_1'))

master <- master %>% left_join(seasonality, by = c('COUNTRY', 'NAME_1'))

# save master dataframe in Github and on M drive
saveRDS(master, './03_output/Senegambia_master.rds') # save
saveRDS(master, 'M:/Hillary/E8_spatial_mixing/Senegambia_master.rds') # save



# E8 map -----------------------------------------------------------------------
# find centroids in E8 admin1s
centroids_E8 <- sf::st_centroid(E8)

# plot
ggplot() +
  geom_sf(data = SSA, fill = "cornsilk2", color = "cornsilk3") +
  geom_sf(data = SADC, fill = "lightgreen", color = "cornsilk3") +
  geom_sf(data = E8, fill = "darkgreen") +
  geom_sf(data = centroids_E8, color = 'red', shape = 4, size = 2) +
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

# plot
ggplot() +
  geom_sf(data = Senegambia_neighbors, fill = "cornsilk2", color = "cornsilk3") +
  geom_sf(data = Senegambia[Senegambia$ID_0 == 'GMB',], fill = "#72B000", color = "cornsilk3") +
  geom_sf(data = Senegambia[Senegambia$ID_0 == 'SEN',], fill = "#72B000", color = "cornsilk3") +
  geom_raster(data = SG_rast, aes(x = x, y = y, alpha = layer), fill = "darkgreen", show.legend = F) +
  geom_sf(data = centroids_combo, aes(color = type), shape = 4, size = 2) +
  theme_bw(base_size = 14) +
  scale_alpha_continuous(range = c(0, 300), breaks = c(0, 0.1, 0.2, 0.5, 1, 300)) +
  scale_color_manual(values = c('red', 'blue')) +
  scale_x_continuous(limits = c(-17.8, -11)) +
  scale_y_continuous(limits = c(12, 17)) +
  labs(x = '', y = '', color = '') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "#daeff8", color = NA))

ggsave('./03_output/Senegambia.pdf', width = 6, height = 6)



# Trip duration ----------------------------------------------------------------
# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-016-1252-3#MOESM2
# 12936_2016_1252_MOESM2_ESM.xlsx Additional file 2, Marshall et al. 2016
library(readxl)
library(splines)

# function to import Marshall et a. 2016 Trips data
import_xlsx <- function(sheet){

  read_excel('./01_data/Marshall_2016_additional_file_2.xlsx', sheet = sheet) %>%
    mutate(Age = ifelse(Age == 'NA', NA, Age),
           Age = is.numeric(Age)) %>%
    dplyr::select(TripID, PersonID, TripType, OrigCountry,
                  DestCountry, Distance, Duration, Purpose, Gender, Age)

}

# import and merge
dat <- map_dfr(.x = c('ML Trips', 'BF Trips', 'ZM Trips', 'TZ Trips'),
               .f = import_xlsx)

# remove NAs (-1s? a little unclear)
dat <- dat %>% mutate(Duration = ifelse(Duration < 0, NA, Duration))
summary(dat$Duration)

# plot by country
ggplot(data = dat,
       aes(x = Distance, y = Duration)) +
  geom_point(alpha = 0.1,  size = 0.5, color = 'darkgrey') +
  geom_smooth(method = 'loess', aes(color = OrigCountry, group = OrigCountry)) +
  geom_smooth(method = 'loess', color = 'black') +
  labs(x = 'Distance (km)',
       y = 'Duration (days)',
       color = 'Origin Country') +
  scale_y_continuous(limits = c(0, 365)) +
  theme_classic()

ggsave('./03_output/Trip_duration.pdf', width = 6, height = 4)


# fit a model to duration data
# linear
model <- glm(Duration ~ Distance, family = gaussian(link = 'identity'), data = dat)
model$coefficients; model$aic

# cubic splines
model <- glm(Duration ~ splines::ns(Distance, 2), family = gaussian(link = 'identity'), data = dat)
model$coefficients; model$aic

# cubic splines
model <- glm(Duration ~ splines::ns(Distance, 6), family = gaussian(link = 'identity'), data = dat)
model$coefficients; model$aic

saveRDS(model, file = './03_output/dist_duration_sk6_model.rds')

# plot by model fit
ggplot(data = dat,
       aes(x = Distance, y = Duration)) +
  geom_point(alpha = 0.1,  size = 0.5, color = 'darkgrey') +
  geom_smooth(formula = y ~ x, method = 'glm', aes(color = 'linear')) +
  geom_smooth(formula = y ~ ns(x, 2), method = 'glm', aes(color = 'spline 2k')) +
  geom_smooth(formula = y ~ ns(x, 6), method = 'glm', aes(color = 'spline 6k')) +
  geom_smooth(method = 'loess', aes(color = 'loess'), size = .5) +
  scale_y_continuous(limits = c(0, 50)) +
  scale_color_discrete(type = c('#00BFC4', 'black', '#7CAE00', '#F8766D')) +
  labs(x = 'Distance (km)',
       y = 'Duration (days)',
       color = 'Fit') +
  theme_classic()

ggsave('./03_output/Trip_duration_fit.pdf', width = 6, height = 4)




# Weighted matrix --------------------------------------------------------------
data <- readRDS('./03_output/Senegambia_master.rds')
model <- readRDS('./03_output/dist_duration_sk6_model.rds')

population <- data$wpop_pop
EIR <- data$EIR

# define parameters [Marshall et al. 2018]
a <- 1.91
log_p <- 4.29
t <- 1.22

# destination population matrix
Nj <- t(matrix(data = c(rep(population, nrow(data))),
               nrow = nrow(data),
               ncol = nrow(data)))

# distance matrix - in km
d <- st_distance(data$w_cent_geometry, data$w_cent_geometry) / 1000

clean_units <- function(x){
  attr(x, "units") <- NULL
  class(x) <- setdiff(class(x), "units")
  x
}

d <- clean_units(d)

# calculate probability of travel to each destination
totals <- (Nj ^ t) * (1 + (d / exp(log_p)))^(-a)

# standardize so that each location sums to 1
Pij <- totals / rowSums(totals)



# add in duration of travel
f <- as_tibble(d) %>% # convert distance matrix to long tibble
  pivot_longer(cols = everything(), values_to = 'Distance', names_to = NULL)

pred <- predict(model, newdata = f) # predict duration of travel with Marshall et al. model fits

pred_dat <- f %>% cbind(pred) %>% # bind predications to data
  mutate(pred = ifelse(Distance == 0, NA, pred)) # remove pred for admin1 with itself

summary(pred_dat$pred) # examine predictions

pred_mat <- matrix(pred_dat$pred, ncol = 20, nrow = 20) # transform predictions back into matrix
pred_mat <- Pij * pred_mat # multiply duration of travel by probability of travel

# replace admin1 with itself relationships with remainder of duration time left to travel
rsum <- 365 - rowSums(pred_mat, na.rm = T)

I <- diag(1, ncol = 20, nrow = 20)
pred_mat[I == 1] <- rsum

# multiply frequency of travel by probability of travel to each destination
# diagonals should always be the majority of the row (staying within own population)
Pij_f <- pred_mat / rowSums(pred_mat)

saveRDS(Pij_f, 'M:/Hillary/E8_spatial_mixing/mixing_W.rds')




