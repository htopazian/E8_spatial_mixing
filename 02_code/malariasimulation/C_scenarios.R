
library(tidyverse)
library(netz) # devtools::install_github('https://github.com/mrc-ide/netz')

data <- readRDS('./03_output/E8_master.rds')


E8 <- c('AGO','BWA','SWZ','MOZ','NAM','ZAF','ZMB','ZWE')
# first, we specify the coutry iso codes and target usages
input <- expand_grid(iso3 = c(c('AGO','BWA','SWZ','MOZ','NAM','ZAF','ZMB','ZWE')),
                     usage = seq(0.1, 0.9, 0.1))

# we can link these to half life, usage rates and distribution frequencies
input <- input %>%
  left_join(get_halflife_data(), by = "iso3") %>%
  left_join(get_usage_rate_data(), by = "iso3") %>%
  mutate(distribution_freq = 365 * 3)

# next we can follow the conversion chain to move between different bed net measures
output <- input %>%
  mutate(access = usage_to_access(usage, usage_rate),
         npc = access_to_crop(access),
         anpcd = crop_to_distribution(npc, distribution_freq = distribution_freq, net_loss_function = net_loss_map, half_life = half_life))
