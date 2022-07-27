
# Specify HPC options ----------------------------------------------------------
library(didehpc)
setwd("M:/Hillary/E8_spatial_mixing/")

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "htopazia")

# transfer the new malariasimulation folder manually to contexts or delete and re-install using conan
# remotes::install_github('mrc-ide/malariasimulation@dev', force=T)
src <- conan::conan_sources(c("https://github.com/mrc-ide/malariasimulation@dev", "https://github.com/mrc-ide/cali"))

ctx <- context::context_save(path = "Q:/contexts",
                             sources = c('M:/Hillary/E8_spatial_mixing/function_PR_EIR_match.R'),
                             packages = c("dplyr", "malariasimulation", "cali"),
                             package_sources = src)

share <- didehpc::path_mapping("malaria", "M:", "//fi--didef3.dide.ic.ac.uk/malaria", "M:")
config <- didehpc::didehpc_config(shares = share,
                                  use_rrq = FALSE,
                                  cores = 1,
                                  cluster = "fi--didemrchnb", # fi--dideclusthn OR fi--didemrchnb
                                  parallel = FALSE)

obj <- didehpc::queue_didehpc(ctx, config = config)


# Set up your job --------------------------------------------------------------
# make all combinations of baseline scenarios
library(tidyverse)
library(malariasimulation)

# read in master Senegambia dataframe
master <- readRDS('M:/Hillary/E8_spatial_mixing/E8_master.rds')

# extra parameters
year <- 365
month <- year/12
human_population <- 10000
sim_length <- 9 * year


# NOTE!!! Need to update species compositions, treatment, etc.

# function to get parameters for baseline scenarios
getparams_calibrate <- function(x){

  p <- master[x,]

  year <- 365
  month <- year/12

  # get starting parameters ----------
  params <- get_parameters(list(
    human_population = human_population,
    model_seasonality = TRUE,
    # rainfall fourier parameters
    # seasonal profiles: c(g0, g[1], g[2], g[3], h[1], h[2], h[3])
    g0 = unlist(p$seasonality)[1],
    g = unlist(p$seasonality)[2:4],
    h = unlist(p$seasonality)[5:7],
    individual_mosquitoes = FALSE))

  # outcome definitions ----------
  params$prevalence_rendering_min_ages = 2 * year
  params$prevalence_rendering_max_ages = 10 * year
  params$clinical_incidence_rendering_min_ages = c(0, 0)*year
  params$clinical_incidence_rendering_max_ages = c(5, 200)*year

  # demography ----------
  flat_demog <- read.table('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing/01_data/Flat_demog.txt') # from mlgts
  ages <- round(flat_demog$V3 * year) # top of age bracket
  deathrates <- flat_demog$V5 / 365 # age-specific death rates

  params <- set_demography(
    params,
    agegroups = ages,
    timesteps = 1,
    deathrates = matrix(deathrates, nrow = 1),
    birthrates = find_birthrates(human_population, ages, deathrates)
  )

  # vectors ----------
  params <- set_species(
    parameters = params,
    species = list(arab_params, fun_params, gamb_params),
    proportions = c(0.25, 0.25, 0.5))

  # proportion of bites taken in bed for each species
  # find values in S.I. of 10.1038/s41467-018-07357-w Table 3
  params$phi_bednets <- c(0.9, 0.9, 0.89) # Hogan et al. 2020
  # proportion of bites taken indoors for each species
  params$phi_indoors <- c(0.96, 0.98, 0.97) # Hogan et al. 2020

  # ITNs ----------
  # find values in S.I. of 10.1038/s41467-018-07357-w
  # or in Table S1.3 of Ellie's 2021 paper
  # same value for all species

  # no resistance
  dn0_1 <- 0.387 # pyr, 0 resistance
  rn_1 <- 0.563 # pyr, 0 resistance
  gamman_1 <- 2.64 # pyr, 0 resistance

  ITNdist <- c(seq(1, (sim_length), 3*year))

  params <- set_bednets(
    parameters = params,
    timesteps = ITNdist,
    coverages = c(rep(p$Mean_ITN_coverage_rate, length(ITNdist))),
    retention = 3 * year,
    dn0 = matrix(c(rep(dn0_1, length(ITNdist))),
                 nrow = length(ITNdist), ncol = 3),
    rn = matrix(c(rep(rn_1, length(ITNdist))),
                nrow = length(ITNdist), ncol = 3),
    rnm = matrix(c(rep(.24, length(ITNdist))),
                   nrow = length(ITNdist), ncol = 3),
    gamman = c(rep(gamman_1 * 365, length(ITNdist)))
    )

  # treatment ----------
  params <- set_drugs(
    parameters = params,
    list(AL_params, SP_AQ_params))

  params$drug_prophylaxis_scale <- c(10.6, 39.34)
  params$drug_prophylaxis_shape <- c(11.3, 3.40)

  params <- set_clinical_treatment(
    parameters = params,
    drug = 1,
    timesteps = c(1),
    coverages = c(0.4)
  )

  parameters <- data.frame(params = c(0), pfpr = c(0))
  parameters$params <- list(params)
  parameters$pfpr <- p$PfPR_rmean
  parameters$COUNTRY <- p$COUNTRY
  parameters$NAME_1 <- p$NAME_1

  return(parameters)

}

parameters <- map_dfr(1:nrow(master), getparams_calibrate)

# save
saveRDS(parameters, 'M:/Hillary/E8_spatial_mixing/calibration_parameters.rds')



# Run tasks -------------------------------------------------------------------
# < submit isolated calibrations ####
# remove ones that have already been run
index <- readRDS('M:/Hillary/E8_spatial_mixing/calibration_parameters.rds') %>%
  mutate(f = paste0('M:/Hillary/E8_spatial_mixing/PR_EIR_match/PR_EIR_admin1_match', '_', COUNTRY, '_', NAME_1, '.rds'),
         index = row_number(),
         exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) %>%
  filter(exist == 0) %>%
  select(index)

t <- obj$enqueue_bulk(index[1:91,], PR_EIR_match)

# read in results
files <- list.files(path = "M:/Hillary/E8_spatial_mixing/PR_EIR_match", pattern = "PR_EIR_admin1_match_", full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))

# concatenate
match <-  do.call("rbind", dat_list) %>% as_tibble()
summary(match$starting_EIR)

# save EIR estimates
saveRDS(match, "M:/Hillary/E8_spatial_mixing/EIRestimates.rds")



# < submit metapopulation calibration ####
t <- obj$enqueue_bulk('W', PR_EIR_match_metapopulation)

