# libraries
library(sf)
library(tidyverse)
library(patchwork)
library(ICDMM)
library(RColorBrewer)
# rebuild ICDMM
# odin::odin_package(getwd())
# click 'install and rebuild'


# HPC set-up -------------------------------------------------------------------
library(didehpc)
setwd('M:/Hillary/E8_spatial_mixing/')

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "htopazia")

# transfer the new odin folder manually to contexts or delete and reinstall using conan

sources <- c('./functions.R')

ctx <- context::context_save(path = "M:/Hillary/contexts",
                             sources = sources,
                             packages = c("dplyr", "odin", "ICDMM", "dde", "tidyr", "nloptr"))

share <- didehpc::path_mapping('malaria', "M:", "//fi--didef3.dide.ic.ac.uk/malaria", "M:")

# template choices: "GeneralNodes", "12Core", "16Core", "12and16Core", "20Core", "24Core" or "32Core"
config <- didehpc::didehpc_config(shares = share,
                                  use_rrq = FALSE,
                                  cores = 10,
                                  cluster = "fi--didemrchnb",
                                  template = "GeneralNodes",
                                  wholenode = TRUE,
                                  parallel = FALSE)

# obj <- didehpc::queue_didehpc(ctx, config = config, provision = "upgrade")
obj <- didehpc::queue_didehpc(ctx, config = config)


# Set up your job --------------------------------------------------------------

master <- readRDS('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing/03_output/Senegambia_master.rds')



# nloptr COBYLA ----------------------------------------------------------------
# list of algorithms: https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
library(nloptr)
source('M:/Hillary/E8_spatial_mixing/functions.R')

# test function - output will be in: 'M:/Hillary/E8_spatial_mixing/optimization_output'
run_nloptr(budget = 500,
           EIR_vector = c(1, 5),
           mixtype = 'I',
           population = c(500, 600))

run_nloptr(combos[1, 1], combos[1, 2], combos[1, 3], combos[1, 4])


# define combinations
budget <- c(0.001 * sum(master$wpop_pop), seq(100000, sum(master$wpop_pop), 100000), sum(master$wpop_pop))
EIR_vector <- list(master$EIR)
mixtype <- c('I', 'W')
population <- list(master$wpop_pop)

combos <- crossing(budget, EIR_vector, mixtype, population)

# submit jobs
obj <- didehpc::queue_didehpc(ctx, config = config)

t <- obj$enqueue_bulk(combos, run_nloptr)
t$status()





# testing simulated annealing --------------------------------------------------
library(GenSA)
source('M:/Hillary/E8_spatial_mixing/functions.R')

# model <- readRDS(model, file = './03_output/dist_duration_sk6_model.rds')
#
# netz <- readRDS("M:/Eradication/rds/netz_median_use_rate.rds")
# conversion <- readRDS("M:/eradication/rds/conversion_usage_pcnets_median_use_rate.rds") %>% as.data.frame() %>%
#   filter(usage<=ITNuse/100)
# conversion$pc_nets_annual <- round(conversion$pc_nets_annual,3)


# objective function:
eval_f0 <- function(itn_cov, budget, EIR_vector, mix, population){

  # print out the ITN usage values which the alogrithm will try next
  print(paste0('Trying itn_cov: ', paste0(itn_cov, collapse = ', ')))

  year <- 365
  n <- length(EIR_vector)

  # define number intervention categories (2 = people with and without nets)
  num_int <- ifelse(sum(itn_cov == 0), 1, 2)

  out <- run_model_metapop(model = "odin_model_metapop",
                           het_brackets = 5,
                           age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,
                                   3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60,
                                   70, 80),
                           init_EIR = EIR_vector,
                           init_ft = 0.4,
                           mix = mix,
                           time = 9 * year,
                           num_int = num_int,
                           ITN_on = 1,
                           ITN_interval = 3 * year,
                           irs_cov = rep(0, n),
                           itn_cov = itn_cov)

  model <- cbind(out$t, out$prev2to10) %>% as_tibble()
  colnames(model) <- c('t', EIR_vector)

  model <- model %>%
    filter(t >= 6*year & t <= 9*year) %>%
    summarize(across(.cols = everything(), .fn = mean)) %>%
    dplyr::select(-t)

  model <- model %>%
    pivot_longer(cols = everything(), names_to = 'EIR', values_to = 'prev2to10') %>%
    cbind(population) %>%
    mutate(cases = prev2to10 * population)

  # set a high penalty when the constraint is not respected
  penalty <- 10

  # minimize overall number of cases 2-10 year olds
  return(sum(model$cases) + penalty * max(sum(itn_cov * population) - budget, 0))

  # use incidence instead of prev

}


# set function to run SA process
solver <- function(budget, EIR_vector, mix, population){

  output <- GenSA(
    par = c(0.5, 0.5),                  # initial values
    fn = eval_f0,                       # function
    lower = rep(0.001, length(EIR_vector)), # lower boundary
    upper = rep(1, length(EIR_vector)), # upper boundary
    control = list(
      temperature = 10
      ), # temperature

    # defining other function inputs
    budget = budget,
    EIR_vector = EIR_vector,
    mix = mix,
    population = population
  )

  result <- tibble(B = budget, eirs = EIR_vector, itn = output$par, b = output$par * population, population) # result

  return(result)


}

# run function
EIR_vector <- c(1, 5)
population <- c(1000, 700)

mix <- matrix(data = c(0.9, 0.1, 0.1, 0.9),
              nrow = 2,
              ncol = 2)

budget <- 300

set.seed(123)
solver(budget, EIR_vector, mix, population)


# check constraint
sum(output$b)

# nb.stop.improvement Integer
