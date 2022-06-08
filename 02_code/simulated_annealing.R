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
                             packages = c("dplyr", "odin", "ICDMM", "dde", "sf", "tidyr"))

share <- didehpc::path_mapping('malaria', "M:", "//fi--didef3.dide.ic.ac.uk/malaria", "M:")

# template choices: "GeneralNodes", "12Core", "16Core", "12and16Core", "20Core", "24Core" or "32Core"
config <- didehpc::didehpc_config(shares = share,
                                  use_rrq = FALSE,
                                  cores = 1,
                                  cluster = "fi--didemrchnb",
                                  template = "32Core",
                                  parallel = FALSE)

# obj <- didehpc::queue_didehpc(ctx, config = config, provision = "upgrade")
obj <- didehpc::queue_didehpc(ctx, config = config)


# Set up your job --------------------------------------------------------------

# define weighted matrix
data <- readRDS('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing/03_output/Senegambia_master.rds')
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

# distance matrix
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

# add in frequency of travel
f <- matrix(diag(c(rep(360/365, 20))), ncol = 20)
I <- diag(1, ncol = 20, nrow = 20)
f[I == 0] <- 5/365

# multiply frequency of travel by probability of travel to each destination
# diagonals should always be the majority of the row (staying within own population)
Pij_f <- (Pij * f) / rowSums(Pij * f)

# distance / population weighted
mixing_W <- (Pij_f)

# save
saveRDS(mixing_W, 'M:/Hillary/E8_spatial_mixing/mixing_W.rds')


# define combinations
mixing <- c('I', 'W')

master <- readRDS('./03_output/Senegambia_master.rds')
EIR <- master$EIR

itn_cov <- seq(0, 1, 0.1)

test <- crossing(EIR, itn_cov)
test



t <- obj$enqueue(odin_mixing(master))
t$status()








# testing COBYLA ---------------------------------------------------------------
# list of algorithms: https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
library(nloptr)

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

  # minimize overall number of cases 2-10 year olds
  return(sum(model$cases))

  # use incidence instead of prev

}

# inequality constraints:
eval_g0 <- function(itn_cov, budget, EIR_vector, mix, population){

  # spending <= budget
  # spending - budget <= 0
  constraint <- sum(itn_cov * population) - budget

  return(constraint)

}

# optimization function:
solver <- function(budget, EIR_vector, mix, population){

  output <- nloptr(
    x0 = c(0.5, 0.5),                # initial values
    eval_f = eval_f0,                # gradient of function
    # eval_g = eval_g0,              # constraint function
    eval_g_ineq = eval_g0,           # inequality constraint function
    # eval_jac_g_ineq = eval_jac_g0, # jacobian of inequality constraint - (derivative of objective function)
    lb = rep(0, length(EIR_vector)), # lower boundary
    ub = rep(1, length(EIR_vector)), # upper boundary
    opts = list("algorithm" = "NLOPT_LN_COBYLA",
                "xtol_rel" = 0.01), # define algorithm

    # defining other function inputs
    budget = budget,
    EIR_vector = EIR_vector,
    mix = mix,
    population = population
  )

  return(output)

  # local or global, eval_f0 function - plot this

}

# run
EIR_vector <- c(1, 5)
population <- c(1000, 700)

mix <- matrix(data = c(0.9, 0.1, 0.1, 0.9),
              nrow = 2,
              ncol = 2)

budget <- 200

set.seed(123)
output <- solver(budget, EIR_vector, mix, population)
output$solution

# check constraint
sum(output$solution * population)



# testing simulated annealing --------------------------------------------------

model <- readRDS(model, file = './03_output/dist_duration_sk6_model.rds')



