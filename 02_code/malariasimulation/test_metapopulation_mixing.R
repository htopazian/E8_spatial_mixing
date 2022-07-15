# testing malariasimulation, metapopulation function

# load libraries
library(dplyr)
library(ggplot2)
library(purrr)
library(data.table)
library(malariaEquilibrium)
library(malariasimulation)
# devtools::install_github('https://github.com/mrc-ide/malariasimulation/tree/feat/metapopulation', force = T)

# set variables
year <- 365
human_population <- 2000
sim_length <- 30 * year
EIR_vector <- c(1, 2, 5)


# get parameters
ms_parameterize <- function(x){ # index of EIR and ITN vector
  
  params <- get_parameters(list(human_population = human_population,
                                model_seasonality = FALSE,
                                individual_mosquitoes = FALSE))
  params <- set_drugs(params, list(AL_params))
  params <- set_clinical_treatment(params, drug = 1, timesteps = 1, coverages = 0.40)

  params <- set_equilibrium(params, init_EIR = EIR_vector[x])
  
  print(paste0('Set params: ', x))
  
  return(params)
  
}

paramslist <- lapply(seq(1, length(EIR_vector), 1), ms_parameterize)


# run model
set.seed(123)

# isolated
output <- run_metapop_simulation(timesteps = sim_length,
                                 parameters = paramslist,
                                 correlations = NULL,
                                 mixing = diag(length(EIR_vector)))

isolate <- data.table::rbindlist(output, idcol = 'EIR') %>%
  mutate(EIR = c(sort(rep(EIR_vector, sim_length))))

# perfectly mixed
output <- run_metapop_simulation(timesteps = sim_length,
                                 parameters = paramslist,
                                 correlations = NULL,
                                 mixing = matrix(rep(.33, 9),
                                                 nrow = 3, ncol = 3))

mixed <- data.table::rbindlist(output, idcol = 'EIR') %>%
  mutate(EIR = c(sort(rep(EIR_vector, sim_length))))

# slightly mixed
output <- run_metapop_simulation(timesteps = sim_length,
                                 parameters = paramslist,
                                 correlations = NULL,
                                 mixing = matrix(c(0.8, 0.1, 0.1, 
                                                   0.1, 0.8, 0.1,
                                                   0.1, 0.1, 0.8),
                                                 nrow = 3, ncol = 3))

middle <- data.table::rbindlist(output, idcol = 'EIR') %>%
  mutate(EIR = c(sort(rep(EIR_vector, sim_length))))


# merge dataframes
mdat <- isolate %>% mutate(model = 'isolated') %>%
  full_join(mixed %>% mutate(model = 'perfectly mixed')) %>%
  full_join(middle %>% mutate(model = 'lightly mixed')) %>%
  mutate(prev2to10 = p_detect_730_3650 / n_730_3650) %>%
  mutate(year = ceiling(timestep/year)) %>%
  group_by(model, EIR, year) %>%
  mutate(prev2to10 = mean(prev2to10, na.rm = T))

ggplot(data = mdat) + 
  geom_line(aes(x = year, y = prev2to10, color = factor(EIR))) +
  facet_wrap(~ model) +
  labs(x = 'time (years)',
       y = 'PfPR 2-10 (month)',
       title = 'Mixing patterns, no ITNs',
       color = 'EIR') +
  theme_classic()

ggsave('C:/Users/htopazia/Desktop/test_mixing.pdf', height = 4, width = 6)

