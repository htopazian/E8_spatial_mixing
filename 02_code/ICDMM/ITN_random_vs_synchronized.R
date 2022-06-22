# Comparing randomly distributed ITNs vs. synchronized distribution between admin1s

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
                                  cores = 4,
                                  cluster = "fi--didemrchnb",
                                  template = "24Core",
                                  # wholenode = TRUE,
                                  parallel = FALSE)

obj <- didehpc::queue_didehpc(ctx, config = config)



# Set up your job
ITN <- c('standard', 'random')
day <- seq(1, 365, 2)
combo <- crossing(ITN, day) %>% filter(!(ITN == 'random' & day > 100))

# submit function to HPC
t <- obj$enqueue_bulk(combo, ICDMM_ITN_test)
t$status()



# Process ----------------------------------------------------------------------
library(data.table)

# read in results
# read in all parameter draw runs and process
files <- list.files(path = "M:/Hillary/E8_spatial_mixing/ITN_timing/", pattern = "ITN_timing*", full.names = TRUE)

dat_list <- lapply(files, function (x) readRDS(x))
output <- rbindlist(dat_list, fill = TRUE, idcol="identifier")

# save output
standard <- output %>% filter(ITN == 'standard')

random <- output %>% filter(ITN == 'random') %>%
  group_by(admin_vector) %>%
  summarize(median = quantile(inc_annual, prob = 0.50),
            q05 = quantile(inc_annual, prob = 0.05),
            q95 = quantile(inc_annual, prob = 0.95))

xrange <- c(1, 365)

random <- crossing(xrange, random)

# plot
ggplot() +
  geom_line(data = random, aes(x = xrange, y = median, color = "random (95% CI)"), size = 1) +
  geom_ribbon(data = random, aes(x = xrange, ymin = q05, ymax = q95), fill = "cornflower blue", color = NA, size = 1, alpha = 0.2) +
  geom_line(data = standard, aes(x = day_seed, y = inc_annual, color = "standardized"), size = 1) +
  scale_y_continuous(limits = c(0, 1.2)) +
  scale_color_manual(values = c("cornflower blue", "blue")) +
  facet_wrap(~ admin_vector) +
  labs(x = 'Day of ITN distribution',
       y = 'Incidence (annual, all ages)',
       color = 'ITN distribution method') +
  theme_classic()

# save
ggsave('./03_output/ICDMM_ITN.pdf', height = 7, width = 15)


# print stats
master <- readRDS('M:/Hillary/E8_spatial_mixing/Senegambia_master.rds')


# find day with lowest resulting incidence
standard %>%
  left_join(master %>% dplyr::select(wpop_pop, NAME_1), by = c('admin_vector' = 'NAME_1')) %>%
  mutate(cases = inc_annual * wpop_pop) %>%
  group_by(day_seed) %>%
  summarize(cases = sum(cases)) %>%
  filter(cases == min(cases)) # day 165, cases 168514


random %>% filter(xrange == 1) %>%
  left_join(master %>% dplyr::select(wpop_pop, NAME_1), by = c('admin_vector' = 'NAME_1')) %>%
  mutate(cases = median * wpop_pop,
         cases05 = q05 * wpop_pop,
         cases95 = q95 * wpop_pop) %>%
  summarize(cases = sum(cases),
            cases05 = sum(cases05),
            cases95 = sum(cases95)) # 205653

205653 / 168514
175617 / 168514
300557 / 168514


