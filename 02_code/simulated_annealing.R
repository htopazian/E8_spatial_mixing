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
                                  cores = 8,
                                  cluster = "fi--didemrchnb",
                                  template = "24Core",
                                  # wholenode = TRUE,
                                  parallel = FALSE)

# obj <- didehpc::queue_didehpc(ctx, config = config, provision = "upgrade")
obj <- didehpc::queue_didehpc(ctx, config = config)


# Set up your job --------------------------------------------------------------

master <- readRDS('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing/03_output/Senegambia_master.rds')



# nloptr COBYLA ----------------------------------------------------------------
# < set up HPC ----
# list of algorithms: https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
library(nloptr)
source('M:/Hillary/E8_spatial_mixing/functions.R')

# test function - output will be in: 'M:/Hillary/E8_spatial_mixing/optimization_output'
# run_nloptr(budget = 500,
#            EIR_vector = c(1, 5),
#            mixtype = 'I',
#            population = c(500, 600))
#
# run_nloptr(combos[1, 1], combos[1, 2], combos[1, 3], combos[1, 4])


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


# < process runs ----
library(data.table)

# list all files run by HPC
files <- list.files(path = "M:/Hillary/E8_spatial_mixing/optimization_output/", pattern = "nloptr_*", full.names = TRUE)

# read in files and combine
dat_list <- lapply(files, function (x) readRDS(x))
dat <- rbindlist(dat_list, fill = TRUE, idcol="file")

# see if constraint is met
dat %>% group_by(B) %>% summarize(b = sum(b)) %>% print(n = 120)

dat2 <- dat %>%
  mutate(b = ifelse(B == 0, 0, b),
         itn = ifelse(B == 0, 0, itn)) %>%
  group_by(B) %>% mutate(t = sum(b)) %>%
  filter(t <= B)



dat2 <- dat2 %>% left_join(master %>% dplyr::select(COUNTRY, PfPR_rmean, EIR), by = c('eirs' = 'EIR'))


# plot
ggplot(data = dat2) +
  geom_line(aes(x = B, y = itn, color = eirs, group = eirs), size = 1, alpha = 0.5) +
  labs(x = '$ USD', y = 'ITN usage', color = 'EIR',
       caption = '1 net = $1 USD, linear fit') +
  # facet_wrap(~COUNTRY, nrow = 1) +
  theme_classic()

ggsave('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing/03_output/nloptr_EIR.pdf', height = 4, width = 6)


dat3 <- dat2 %>% filter(B == 1000000) %>%
  mutate(start = B / sum(population)) %>%
  mutate(diff = itn - start)

head(dat3)

summary(dat3$itn); summary(dat3$diff)



countries <- readRDS('./03_output/countries.rds')


Senegambia <- c('SEN','GMB')

Senegambia_neighbors <- c('GNB','GIN','MRT','MLI')

Senegambia <- countries %>%
  filter(ID_0 %in% Senegambia)

Senegambia_neighbors <- countries %>%
  filter(ID_0 %in% Senegambia_neighbors)

ggplot() +
  geom_sf(data = Senegambia_neighbors, fill = "cornsilk2", color = "cornsilk3") +
  geom_sf(data = Senegambia[Senegambia$ID_0 == 'GMB',], fill = "#72B000", color = "cornsilk3") +
  geom_sf(data = Senegambia[Senegambia$ID_0 == 'SEN',], fill = "#72B000", color = "cornsilk3") +
  geom_sf(data = st_as_sf(dat3), aes(fill = diff)) +
  theme_bw(base_size = 14) +
  # https://colorbrewer2.org/#type=diverging&scheme=RdBu&n=3
  scale_fill_gradient2(low = '#67a9cf', mid = "white", high = '#ef8a62', midpoint = 0) +
  scale_x_continuous(limits = c(-17.8, -11)) +
  scale_y_continuous(limits = c(12, 17)) +
  labs(x = '', y = '', fill = 'change in ITN use \nafter optimization',
       caption = '1 M USD, weighted mixing, vs. a constant usage of 0.462') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "#daeff8", color = NA),
        plot.caption.position = "panel",
        plot.caption = element_text(hjust = 0))

ggsave('./03_output/Senegambia_optimization_diff.pdf', width = 6, height = 6)





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
    par = budget * (population / sum(population)) / population, # initial values
    fn = eval_f0,                       # function
    lower = rep(0.001, length(EIR_vector)), # lower boundary
    upper = rep(1, length(EIR_vector)), # upper boundary
    control = list(
      temperature = .1,
      simple.function = TRUE,
      acceptance.param = 2,
      nb.stop.improvement = 2
      # max.call = 2
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

# NETZ ----
library(netz)

devtools::install_github('https://github.com/mrc-ide/netz')

# specify the coutry iso codes and target usages
input <- expand_grid(iso3 = c("SEN", "GMB"),
                     usage = seq(0, 1, 0.01))

# link these to half life, usage rates and distribution frequencies
input <- input %>%
  left_join(get_halflife_data(), by = "iso3") %>%
  left_join(get_usage_rate_data(), by = "iso3") %>%
  mutate(distribution_freq = 365 * 3)

# follow the conversion chain to move between different bed net measures
output <- input %>%
  mutate(access = usage_to_access(usage, usage_rate),
         npc = access_to_crop(access, type = "loess_extrapolate"),
         anpcd = crop_to_distribution(npc, distribution_freq = distribution_freq, net_loss_function = net_loss_map, half_life = half_life))

# Plot curve
ggplot() +
  geom_line(data = output, aes(x = anpcd, y = usage, color = iso3), size = 1.5) +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(x = "Annual nets distributed per capita", y = "Usage", colour = "") +
  theme_classic()


ggsave('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing/03_output/netz.pdf', height = 4, width = 4)

