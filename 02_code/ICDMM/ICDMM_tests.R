# Testing ICDMM
# Running 20 admin units on HPC as isolated models (no-mixing)
# Running 20 admin units individually using run_model
# Comparing results


# Libraries
library(tidyverse)
library(patchwork)
library(ICDMM)
# rebuild ICDMM
# odin::odin_package(getwd())
# click 'install and rebuild'


# Run admins together ----------------------------------------------------------
# HPC set-up
library(didehpc)
setwd('M:/Hillary/E8_spatial_mixing/')

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "htopazia")

# transfer the new odin folder manually to contexts or delete and reinstall using conan
# C:\Users\htopazia\Documents\R\win-library\4.1

sources <- c('./functions.R')

ctx <- context::context_save(path = "M:/Hillary/contexts",
                             sources = sources,
                             packages = c("dplyr", "odin", "ICDMM", "dde", "tidyr", "nloptr"))

share <- didehpc::path_mapping('malaria', "M:", "//fi--didef3.dide.ic.ac.uk/malaria", "M:")

# template choices: "GeneralNodes", "12Core", "16Core", "12and16Core", "20Core", "24Core" or "32Core"
config <- didehpc::didehpc_config(shares = share,
                                  use_rrq = FALSE,
                                  cores = 8,
                                  cluster = "fi--didemrchnb",
                                  template = "24Core",
                                  # wholenode = TRUE,
                                  parallel = FALSE)

obj <- didehpc::queue_didehpc(ctx, config = config)



# Set up your job
# read in master Senegambia dataset
master <- readRDS('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing/03_output/Senegambia_master.rds')

# submit function to HPC
t <- obj$enqueue(ICDMM_admin_HPC(master))
t$status()

# read in results
model1 <- readRDS('./ICDMM_tests.rds')



# Run admins individually ------------------------------------------------------
setwd('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing')

age_vector <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5,
                7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
year <- 365

country_vector <- master$COUNTRY
admin_vector <- master$NAME_1
EIR_vector <- master$EIR

n <- length(EIR_vector)

ITN_vector <- seq(year, year*n, year)

combinations <- as_tibble(cbind(EIR_vector, country_vector, admin_vector, ITN_vector)) %>%
  mutate(EIR_vector = as.numeric(EIR_vector),
         ITN_vector = as.numeric(ITN_vector))


ICDMM_admin <- function(x){

  combo <- combinations[x, ]

  init_EIR <- combo$EIR_vector
  country <- combo$country_vector
  admin2 <- combo$admin_vector
  ITN_IRS_on <- combo$ITN_vector

  out <- run_model(model = "odin_model",
                   het_brackets = 5,
                   age = age_vector,
                   init_EIR = init_EIR,
                   init_ft = 0.4,
                   country = country,
                   admin2 = admin2,
                   time = 25 * year,
                   num_int = 2,
                   ITN_IRS_on = ITN_IRS_on,
                   ITN_interval = 3 * year,
                   itn_cov= 0.5,
                   irs_cov = 0)


  model <- cbind(out$t, out$prev2to10) %>% as_tibble()
  colnames(model) <- c('t', admin2)
  model <- model %>%
    pivot_longer(cols = 2, names_to = 'admin_vector', values_to = 'prev2to10')

  return(model)

}

index <- seq(1, nrow(combinations), 1)

model2 <- map_dfr(index, ICDMM_admin)

saveRDS(model2, 'M:/Hillary/E8_spatial_mixing/ICDMM_tests2.rds')



# Compare ----------------------------------------------------------------------
# merge

modeldat <- model1 %>% mutate(model = 'aggregated run') %>%
  full_join(model2 %>% mutate(model = 'isolated runs') )


segments <- tibble(admin_vector, ITN_vector)

ggplot(data = modeldat) +
  geom_vline(data = segments, aes(xintercept = ITN_vector/365), color = 'lightgrey', lty = 2) +
  geom_line(aes(x = t/365, y = prev2to10, color = model, size = model)) +
  scale_y_continuous(limits = c(0, 0.12)) +
  scale_color_manual(values = c("cornflower blue", "blue")) +
  scale_size_manual(values = c(1.5, .2), guide = 'none') +
  facet_wrap(~ admin_vector) +
  labs(x = 'time (years)',
       y = 'PfPR 2-10',
       color = '') +
  theme_classic()


ggsave('./03_output/ICDMM_test_20.pdf', height = 7, width = 15)
