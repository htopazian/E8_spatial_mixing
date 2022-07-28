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
# E8
E8 <- c('AGO','BWA','SWZ','MOZ','NAM','ZAF','ZMB','ZWE')

# SADC
SADC <- c('COM','COD','LSO','MDG','MWI','MUS','SYC','TZA')

# other countries in SSA
SSA <- c('SEN','GMB','GNB','GIN','MRT','MLI','KEN','UGA','COG','SSD','ETH','SOM','CAF','CMR','GAB','GNQ','RWA','BDI','NGA','TCD')


# saves copy of gadm SpatVector as .rds in data folder
import_gadm <- function(ISO, level){
  geodata::gadm(country = ISO, level = level, path = './01_data/gadm', version="4.0")
}

map2(E8, 1, import_gadm)    # admin1
map2(SADC, 0, import_gadm)  # admin1
map2(SSA, 0, import_gadm)   # admin0


# create a list of all countries
# E8 and SADC admin1s
admin1 <- lapply(c(E8), function(x){
    list.files(path = "./01_data/gadm",
               pattern = paste0("*", x , "_1_pk.rds"), full.names = TRUE)
  }) %>% unlist

# SSA countries admin0s
admin0 <- lapply(c(SSA, SADC), function(x){
  list.files(path = "./01_data/gadm",
             pattern = paste0("*", x , "_0_pk.rds"), full.names = TRUE)
}) %>% unlist

files <- c(admin1, admin0)

# unpack gadms
unpack_gadm <- function(file){

  object <- readRDS(file) # read in object
  object <- terra::vect(object) # unpack SpatVector
  st_as_sf(object) # transform to sf object

}

countries <- map_dfr(files, unpack_gadm) # loop over each country
st_crs(countries) # view CRS

saveRDS(countries, './03_output/countries.rds') # save
# countries <- readRDS('./03_output/countries.rds')


# subset countries to E8
E8_sf <- countries %>% filter(ID_0 %in% E8)


# Worldpop ---------------------------------------------------------------------
# https://github.com/wpgp/wpgpDownloadR

# download raster file and .csv of population counts
lapply(E8, function(x){

  wpgpGetCountryDataset(ISO3 = x, covariate = 'ppp_2020', destDir = './01_data/worldpop')
  wpgpGetPOPTable(ISO3 = x, year = 2020, destDir = './01_data/worldpop')

  print(paste0('Read in population for country ', x))

})

# import rasters and merge
E8_rast <- raster::merge(raster('./01_data/worldpop/ago_ppp_2020.tif'),
                         raster('./01_data/worldpop/bwa_ppp_2020.tif'),
                         raster('./01_data/worldpop/swz_ppp_2020.tif'),
                         raster('./01_data/worldpop/moz_ppp_2020.tif'),
                         raster('./01_data/worldpop/nam_ppp_2020.tif'),
                         raster('./01_data/worldpop/zaf_ppp_2020.tif'),
                         raster('./01_data/worldpop/zmb_ppp_2020.tif'),
                         raster('./01_data/worldpop/zwe_ppp_2020.tif'))
# check res
res(E8_rast)

saveRDS(E8_rast, './03_output/E8_rast.rds') # save
# E8_rast <- readRDS('./03_output/E8_rast.rds')


# downsize resolution for ease of plotting
E8_rast <- raster::aggregate(E8_rast, fact = 100); res(E8_rast)

# extract population within each admin1
pop <- raster::extract(E8_rast,
                       E8_sf,
                       exact = T)

# sum for each polygon
pop_2020 <- unlist(lapply(pop, function(x) if (!is.null(x)) sum(x, na.rm = TRUE) else NA)) %>%
  as_tibble() %>%
  rename(wpop_pop = value)

pop_2020 <- tibble(E8_sf$ID_0, E8_sf$NAME_1, pop_2020)

colnames(pop_2020) <- c('ID_0', 'NAME_1', 'wpop_pop')


# change to dataframe to use with ggplot()
E8_rast_reduced <- data.frame(rasterToPoints(E8_rast, spatial = T))
saveRDS(E8_rast_reduced, './03_output/E8_rast_reduced.rds') # save


# MAP ITN ----------------------------------------------------------------------
ITN_2020 <- read_csv('./01_data/map/00_ITN_table_Africa_admin1_2000-2020.csv') %>%
  filter(ISO %in% E8 & Year == 2020) %>%
  dplyr::select(ISO, Name_1, Mean_ITN_coverage_rate) %>%
  mutate(Name_1 = case_when(ISO == 'AGO' & Name_1 == 'Bie' ~ 'Bié',
                            ISO == 'AGO' & Name_1 == 'Huila' ~ 'Huíla',
                            ISO == 'AGO' & Name_1 == 'Uige' ~ 'Uíge',
                            ISO == 'MOZ' & Name_1 == 'Niassa' ~ 'Nassa',
                            ISO == 'BWA' & Name_1 == 'North East' ~ 'North-East',
                            ISO == 'BWA' & Name_1 == 'North West' ~ 'North-West',
                            ISO == 'BWA' & Name_1 == 'South East' ~ 'South-East',
                            TRUE ~ Name_1))


# MAP PfPR ---------------------------------------------------------------------
# https://malariaatlas.org/malaria-burden-data-download/
# Annual mean of PF Parasite Rate
pfpr <- read_csv('./01_data/map/00_PfPR_table_Global_admin1_2000-2019.csv') %>%
  filter(Year == 2019 & ISO %in% E8) %>%
  dplyr::select(ISO, Name_0, Name_1, Pop, PfPR_rmean) %>% rename(map_pop = Pop) %>%
  mutate(Name_1 = case_when(ISO == 'AGO' & Name_1 == 'Bie' ~ 'Bié',
                            ISO == 'AGO' & Name_1 == 'Huila' ~ 'Huíla',
                            ISO == 'AGO' & Name_1 == 'Uige' ~ 'Uíge',
                            ISO == 'MOZ' & Name_1 == 'Niassa' ~ 'Nassa',
                            ISO == 'BWA' & Name_1 == 'North East' ~ 'North-East',
                            ISO == 'BWA' & Name_1 == 'North West' ~ 'North-West',
                            ISO == 'BWA' & Name_1 == 'South East' ~ 'South-East',
                            TRUE ~ Name_1))


# Centroids --------------------------------------------------------------------
# find centroids in E8 admin1s
centroids <- sf::st_centroid(E8_sf) %>%
  dplyr::select(ID_0, NAME_1, geometry) %>%
  rename(cent_geometry = geometry) %>% as_tibble()

saveRDS(centroids, './03_output/E8_centroids.rds')

# weighted polygon centroids
# convert polygons to raster layer
z <- rasterize(E8_sf, E8_rast)
# compute weighted x and y coordinates within each rasterized region
xx <- zonal(init(E8_rast, v = "x") * E8_rast, z) / zonal(E8_rast, z) # default = mean
yy <- zonal(init(E8_rast, v = "y") * E8_rast, z) / zonal(E8_rast, z) # default = mean

# combine results in a matrix
weighted_centroids <- tibble(xx[, 2], yy[, 2])

weighted_centroids <- st_as_sf(weighted_centroids, coords = c('xx[, 2]', 'yy[, 2]'))

weighted_centroids <- tibble(E8_sf$ID_0, E8_sf$NAME_1, weighted_centroids)
colnames(weighted_centroids) <- c('ID_0', 'NAME_1', 'w_cent_geometry')

saveRDS(weighted_centroids, './03_output/E8_weighted_centroids.rds')



# Seasonality ------------------------------------------------------------------
# specify time-range
start_date <- "2020-01-01"
end_date <- "2020-12-31"

# function to pull and process rainfall data using the umbrella package
process_rain <- function(x){ # x = index of admin1

  E8_subset <- E8_sf[x,]

  message(paste("Pulling rainfall for", E8_subset$ID_0, E8_subset$NAME_1, sep = " "))

  # extract
  daily_rain_raw <- umbrella::pull_daily_rainfall(
    sf = E8_subset,
    start_date = start_date,
    end_date = end_date)

  # process the raw data
  daily_rain  <- daily_rain_raw %>%
    # convert to long and format
    pivot_longer(-c(ID_0:ISO_1),
                 names_to = "date",
                 values_to = "rainfall",
                 names_prefix = "X",
                 names_pattern = "(.*)_precipitation") %>%
    mutate(date = as.Date(as.character(readr::parse_number(.data$date)), format = "%Y%m%d"),
           year = lubridate::year(.data$date),
           day_of_year = lubridate::yday(.data$date),
           t = lubridate::yday(.data$date) / 365,
           rainfall = as.numeric(rainfall)) %>%
    # replace missing data with 0
    replace_na(replace = list(rainfall = 0)) %>%
    # remove any leap year addtional days
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

  seasonality <- cbind(rain$COUNTRY, rain$NAME_1, seasonality)
  colnames(seasonality) <- c('COUNTRY', 'NAME_1', 'seasonality')

  return(seasonality)

}

# Run function for one admin1 unit at a time.
# Running entire countries at a time requires too much computational power.
seasonality <- map_dfr(c(1:nrow(E8_sf)),
                       process_rain)

saveRDS(seasonality, './03_output/E8_seasonality.rds')



# Save master dataset ----------------------------------------------------------
# join all created datasets
master <- E8_sf %>%
  left_join(pop_2020) %>%
  left_join(ITN_2020, by = c('ID_0' = 'ISO', 'NAME_1' = 'Name_1')) %>%
  left_join(pfpr, by = c('ID_0' = 'ISO', 'NAME_1' = 'Name_1')) %>%
  left_join(centroids) %>%
  left_join(weighted_centroids) %>%
  left_join(seasonality, by = c('COUNTRY', 'NAME_1')) %>%
  dplyr::select(ID_0, COUNTRY, NAME_1, map_pop, wpop_pop, Mean_ITN_coverage_rate,
                PfPR_rmean, seasonality, cent_geometry, w_cent_geometry, geometry)

master %>% filter(is.na(Mean_ITN_coverage_rate)) # 6 missing, all in BWA
master %>% filter(is.na(PfPR_rmean)) # 1 missing in MOZ

# correct missings
master <- master %>%
  mutate(PfPR_rmean = ifelse(NAME_1 == 'Maputo City',
                             master[master$NAME_1 == 'Maputo', ]$PfPR_rmean,
                             PfPR_rmean),
         Mean_ITN_coverage_rate = ifelse(is.na(Mean_ITN_coverage_rate),
                                         0,
                                         Mean_ITN_coverage_rate))

# save master dataframe in Github and on M drive
saveRDS(master, './03_output/E8_master.rds')
saveRDS(master, 'M:/Hillary/E8_spatial_mixing/E8_master.rds')


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
  geom_smooth(method = 'loess', aes(fill = 'Overall fit'), color = 'black') +
  scale_fill_manual(values = c(NA)) +
  labs(x = 'Distance (km)',
       y = 'Duration (days)',
       color = '',
       fill = '') +
  coord_cartesian(ylim = c(0, 365)) +
  theme_classic() +
  theme(text = element_text(size = 14))

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

# loess
model <- loess(Duration ~ Distance, data = dat)
model$coefficients; model$aic

saveRDS(model, file = './03_output/dist_duration_loess_model.rds')

# plot by model fit
ggplot(data = dat,
       aes(x = Distance, y = Duration)) +
  geom_point(alpha = 0.1,  size = 0.5, color = 'darkgrey') +
  geom_smooth(formula = y ~ x, method = 'glm', aes(color = 'linear')) +
  geom_smooth(formula = y ~ ns(x, 2), method = 'glm', aes(color = 'spline 2k')) +
  geom_smooth(formula = y ~ ns(x, 6), method = 'glm', aes(color = 'spline 6k')) +
  geom_smooth(method = 'loess', aes(color = 'loess')) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_color_discrete(type = c('#00BFC4', 'black', '#7CAE00', '#F8766D')) +
  labs(x = 'Distance (km)',
       y = 'Duration (days)',
       color = 'Fit') +
  theme_classic()

ggsave('./03_output/Trip_duration_fit.pdf', width = 6, height = 4)



# Weighted matrix --------------------------------------------------------------
data <- readRDS('./03_output/E8_master.rds')
# model <- readRDS('./03_output/dist_duration_loess_model.rds')
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
  mutate(pred = ifelse(Distance == 0, NA, pred), # remove pred for admin1 with itself
         pred = ifelse(pred <0, 0, pred))        # turn negative durations to 0

summary(pred_dat$pred) # examine predictions

pred_mat <- matrix(pred_dat$pred, ncol = nrow(data), nrow = nrow(data)) # transform predictions back into matrix
pred_mat <- Pij * pred_mat # multiply duration of travel by probability of travel

# replace admin1 with itself relationships with remainder of duration time left to travel
rsum <- 365 - rowSums(pred_mat, na.rm = T)

I <- diag(1, ncol = nrow(data), nrow = nrow(data))
pred_mat[I == 1] <- rsum


# multiply frequency of travel by probability of travel to each destination
# diagonals should always be the majority of the row (staying within own population)
Pij_f <- pred_mat / rowSums(pred_mat)

saveRDS(Pij_f, 'M:/Hillary/E8_spatial_mixing/E8_mixing_W.rds')


# E8 map -----------------------------------------------------------------------
countries <- readRDS('./03_output/countries.rds')
weighted_centroids <- readRDS('./03_output/E8_weighted_centroids.rds') %>% st_as_sf(crs = crs(countries))
E8_rast_reduced <- readRDS('./03_output/E8_rast_reduced.rds')

# base plot
baseplot <- ggplot() +
  geom_sf(data = countries %>% filter(ID_0 %in% SSA), fill = "cornsilk2", color = "cornsilk3", size = 0.01) +
  geom_sf(data = countries %>% filter(ID_0 %in% SADC), aes(fill = 'SADC'), color = "cornsilk3", size = 0.01, alpha = 0.6) +
  geom_sf(data = E8_sf, aes(fill = 'E8'), color = "cornsilk3", alpha = 0.6, size = 0.01) +
  scale_fill_manual(values = c("#55C667FF", "#95D840FF")) +
  theme_bw(base_size = 14) +
  coord_sf(xlim = c(10, 51.5), ylim = c(-35, 6)) +
  labs(x = '', y = '', color = '', fill = '') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "#daeff8", color = NA))

# plot population
# choose points for population plotting (alpha breaks)
summary(E8_rast_reduced$layer)
# create point for legend
x <- c(1); y <- c(1); point <- tibble(x, y)

A <- baseplot +
  geom_tile(data = E8_rast_reduced, aes(x = x, y = y, alpha = layer), fill = "darkgreen", show.legend = F) +
  geom_point(data = point, aes(x = x, y = y, color = "population"), size = 0.5) +
  scale_color_manual(values = c("darkgreen")) +
  scale_alpha_continuous(range = c(0, 420), breaks = c(0, 0.0038, 0.0038, 0.1216, 420))

# plot weighted centroids
B <- baseplot +
  geom_sf(data = weighted_centroids, aes(color = 'weighted centroids'), shape = 4, size = 1) +
  scale_color_manual(values = c('blue')) +
  coord_sf(xlim = c(10, 51.5), ylim = c(-35, 6))


# plot flow between centroids
data <- readRDS('./03_output/E8_master.rds')
matrix <- readRDS('M:/Hillary/E8_spatial_mixing/E8_mixing_W.rds')

test <- as.data.frame(matrix)
test <- test %>% pivot_longer(cols = everything(), names_to = 'destination', values_to = 'flow')

test$country_origin <- unlist(lapply(master$ID_0, function(x){rep(x, nrow(master))}))
test$admin1_origin <- unlist(lapply(master$NAME_1, function(x){rep(x, nrow(master))}))
test$country_destination <- rep(master$ID_0, nrow(master))
test$admin1_destination <- rep(master$NAME_1, nrow(master))

test2 <- test %>%
  left_join(weighted_centroids, by = c('country_origin' = 'ID_0', 'admin1_origin' = 'NAME_1')) %>%
  rename(centroid_origin = w_cent_geometry) %>%
  left_join(weighted_centroids, by = c('country_destination' = 'ID_0', 'admin1_destination' = 'NAME_1')) %>%
  rename(centroid_destination = w_cent_geometry) %>%
  mutate(long_origin = unlist(map(centroid_origin, 1)),
         lat_origin = unlist(map(centroid_origin, 2)),
         long_destination = unlist(map(centroid_destination, 1)),
         lat_destination = unlist(map(centroid_destination, 2))) %>%
  mutate(flow_f = case_when(flow == 0 ~ 'none',
                            flow > 0 & flow <= 0.01 ~ 'less movement',
                            flow > 0.01 & flow <= 0.1 ~ 'more movement',
                            flow > 0.1 ~ 'max movement'),
         flow_f = factor(flow_f)) %>%
  filter(flow > 0 & flow < 0.9)

summary(test2$flow)

C <- ggplot() +
  geom_sf(data = countries %>% filter(ID_0 %in% c(SSA, SADC, E8)), fill = "cornsilk2", color = "cornsilk3", size = 0.01) +
  geom_sf(data = st_union(countries %>% filter(ID_0 %in% SADC)), color = "#95D840FF", fill = NA, size = 0.5) +
  geom_sf(data = st_union(countries %>% filter(ID_0 %in% E8)), color = "#55C667FF", fill = NA, size = 0.5) +
  geom_segment(data = test2, aes(x = long_origin, y = lat_origin, xend = long_destination, yend = lat_destination, alpha = flow_f), size = 0.1, color =  "blue") +
  scale_alpha_manual(values = c(0.05, 0.1), name = '', labels = c('less movement', 'more movement')) +
  theme_bw(base_size = 14) +
  coord_sf(xlim = c(10, 51.5), ylim = c(-35, 6)) +
  labs(x = '', y = '', color = '', fill = '', alpha = '') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "#daeff8", color = NA))


# consolidate plots
mapplot <- A + B + C + plot_layout(guides = "collect", nrow = 1) +
  plot_annotation(tag_levels = 'A')

ggsave('./03_output/E8_map.pdf', plot = mapplot, width = 12, height = 6)
