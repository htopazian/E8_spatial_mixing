# testing malariasimulation, metapopulation function
library(dplyr)
library(ggplot2)
library(data.table)
library(malariaEquilibrium)
library(malariasimulation)

# devtools::install_github('https://github.com/mrc-ide/malariasimulation/tree/feat/metapopulation', force = T)


year <- 365
human_population <- 10000
sim_length <- 10 * year
EIR_vector <- c(1, 2, 5)
ITN_vector <- list(c(5, 6, 7, 8)*year, c(6, 7, 8)*year, c(7, 8)*year)

# seasonalities
highly_seasonal <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
seasonal <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
perennial <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))
seasonality <- c(highly_seasonal, seasonal, perennial)


set.seed(123)

# run individual EIRs
ms_isolated <- function(x){ # index of EIR vector

  params <- get_parameters(list(human_population = human_population,
                                # model_seasonality = TRUE,
                                # # rainfall fourier parameters
                                # g0 = unlist(seasonality[x])[1],
                                # g = unlist(seasonality[x])[2:4],
                                # h = unlist(seasonality[x])[5:7],
                                individual_mosquitoes = FALSE))
  params <- set_drugs(params, list(AL_params))
  params <- set_clinical_treatment(params, drug = 1, timesteps = 1, coverages = 0.40)

  n <- length(ITN_vector[[x]]) # number of ITN distributions

  params <- set_bednets(
    params,
    timesteps = ITN_vector[[x]],
    coverages = rep(0.5, n),
    retention = 3 * year,
    dn0 = matrix(c(.533, .45), nrow = n, ncol = 1),
    rn = matrix(c(.56, .5), nrow = n, ncol = 1),
    rnm = matrix(c(.24, .24), nrow = n, ncol = 1),
    gamman = rep(2.64 * year, n))

  params <- set_equilibrium(params, init_EIR = EIR_vector[x])
  output <- run_simulation(timesteps = sim_length, params)
  out <- output %>% mutate(EIR = EIR_vector[x])

}

model1 <- map_dfr(seq(1, length(EIR_vector), 1), ms_isolated)


# run metapopulation model
ms_mixed <- function(x){

  params <- get_parameters(list(human_population = human_population,
                                # model_seasonality = TRUE,
                                # # rainfall fourier parameters
                                # g0 = unlist(seasonality[x])[1],
                                # g = unlist(seasonality[x])[2:4],
                                # h = unlist(seasonality[x])[5:7],
                                individual_mosquitoes = FALSE))
  params <- set_drugs(params, list(AL_params))
  params <- set_clinical_treatment(params, drug = 1, timesteps = 1, coverages = 0.40)

  n <- length(ITN_vector[[x]]) # number of ITN distributions

  params <- set_bednets(
    params,
    timesteps = ITN_vector[[x]],
    coverages = rep(0.5, n),
    retention = 3 * year,
    dn0 = matrix(c(0.387), nrow = n, ncol = 1),
    rn = matrix(c(0.563), nrow = n, ncol = 1),
    rnm = matrix(c(0.24), nrow = n, ncol = 1),
    gamman = rep(2.64 * year, n))

  params <- set_equilibrium(params, init_EIR = EIR_vector[x])

}

paramslist <- lapply(seq(1, length(EIR_vector), 1), ms_mixed)

output <- run_metapop_simulation(timesteps = sim_length,
                                 parameters = paramslist,
                                 correlations = NULL,
                                 mixing = diag(length(EIR_vector)))


model2 <- data.table::rbindlist(output, idcol = 'EIR') %>%
  mutate(EIR = c(sort(rep(EIR_vector, sim_length))))


mdat <- model1 %>% mutate(model = 'isolated') %>%
  full_join(model2 %>% mutate(model = 'mixed')) %>%
  mutate(prev2to10 = p_detect_730_3650 / n_730_3650) %>%
  mutate(month = ceiling(timestep/(year/12))) %>%
  group_by(model, EIR, month) %>%
  mutate(prev2to10 = mean(prev2to10, na.rm = T))

table(mdat$EIR, mdat$model)


# plot
num_dist <- unlist(lapply(ITN_vector, length))
years_dist <- unlist(ITN_vector) / year

EIR <- unlist(lapply(seq(1, length(EIR_vector), 1), function(x){rep(EIR_vector[x], num_dist[x])}))
segments <- tibble(EIR, years_dist)

ggplot(data = mdat) +
  geom_vline(data = segments, aes(xintercept = years_dist), color = 'lightgrey', lty = 2) +
  geom_line(aes(x = month / 12, y = prev2to10, color = model, size = model)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  scale_color_manual(values = c("cornflower blue", "blue")) +
  scale_size_manual(values = c(2, 1), guide = 'none') +
  facet_wrap(~ EIR, labeller = label_both) +
  labs(x = 'time (years)',
       y = 'PfPR 2-10 (month)',
       title = 'run_simulation vs. run_metapop_simulation',
       # subtitle = 'three seasonalities, three ITN timings',
       subtitle = 'three ITN timings',
       color = '',
       caption = paste0('population = ', human_population)) +
  theme_classic()


ggsave('C:/Users/htopazia/Desktop/test_ITN_seasonality.pdf', height = 4, width = 6)
ggsave('C:/Users/htopazia/Desktop/test_ITN.pdf', height = 4, width = 6)
