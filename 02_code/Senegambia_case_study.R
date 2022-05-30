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
# read in datasets from 'maps.R'
master <- readRDS('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing/03_output/Senegambia_master.rds')

t <- obj$enqueue(odin_mixing(master))
t$status()


# Process ----------------------------------------------------------------------
setwd('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing')

output <- readRDS('M:/Hillary/E8_spatial_mixing/output.rds')

# create PfPR plotting function
mix_model_plot <- function(mix_name){

  ggplot(data = output %>% filter(name == mix_name)) +
    geom_line(aes(x = year, y = prev2to10, group = EIR, color = as.numeric(EIR)), show.legend = T) +
    scale_y_continuous(limits = c(0,0.1)) +
    labs(x = 'time (years)',
         y = 'PfPR 2-10',
         color = 'EIR',
         title = mix_name) +
    theme_classic()

}

# run three scenarios
A <- mix_model_plot('Isolated')
B <- mix_model_plot('Perfectly mixed')
C <- mix_model_plot('Weighted')

A + B + C + plot_layout(guides = "collect", nrow=1) + plot_annotation(tag_levels = 'A')

# save
ggsave('./03_output/Senegambia_mixing_PfPR.pdf', width = 8, height = 4)


# read in datasets from 'maps.R'
master <- readRDS('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing/03_output/Senegambia_master.rds')


countries <- readRDS('./03_output/countries.rds')
Senegambia <- c('SEN','GMB')
Senegambia_neighbors <- c('GNB','GIN','MRT','MLI')

Senegambia <- countries %>%
  filter(ID_0 %in% Senegambia)

Senegambia_neighbors <- countries %>%
  filter(ID_0 %in% Senegambia_neighbors)


# plot isolated PfPR prevalence by admin1
pfpr_diff <- output %>% filter(year == 10 & name %in% c('Isolated')) %>%
  group_by(EIR, name) %>%
  summarize(prev2to10 = mean(prev2to10)) %>%
  dplyr::select(-name) %>%
  mutate(EIR = round(as.numeric(EIR),6))

plot_data <- master %>%
  mutate(EIR = round(as.numeric(EIR),6)) %>%
  left_join(pfpr_diff, by = 'EIR')

ggplot() +
  geom_sf(data = Senegambia_neighbors, fill = "cornsilk2", color = "cornsilk3") +
  geom_sf(data = Senegambia[Senegambia$ID_0 == 'GMB',], fill = "#72B000", color = "cornsilk3") +
  geom_sf(data = Senegambia[Senegambia$ID_0 == 'SEN',], fill = "#72B000", color = "cornsilk3") +
  geom_sf(data = plot_data, aes(fill = prev2to10)) +
  theme_bw(base_size = 14) +
  scale_alpha_continuous(range = c(0, 300), breaks = c(0, 0.1, 0.2, 0.5, 1, 300)) +
  scale_fill_distiller(palette = 'Greens', direction = 1) +
  scale_x_continuous(limits = c(-17.8, -11)) +
  scale_y_continuous(limits = c(12, 17)) +
  labs(x = '', y = '', fill = 'PfPR no mixing') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "#daeff8", color = NA))

ggsave('./03_output/Senegambia_mixing_I.pdf', width = 6, height = 6)


# plot difference between isolated and weighted PfPR prevalence by admin1
pfpr_diff <- output %>% filter(year == 10 & name %in% c('Isolated', 'Weighted')) %>%
  group_by(EIR, name) %>%
  summarize(prev2to10 = mean(prev2to10)) %>%
  arrange(EIR, name) %>%
  mutate(diff = prev2to10 - lag(prev2to10)) %>%
  filter(!is.na(diff)) %>%
  dplyr::select(-name) %>%
  mutate(EIR = round(as.numeric(EIR),6))

plot_data <- master %>%
  mutate(EIR = round(as.numeric(EIR),6)) %>%
  left_join(pfpr_diff, by = 'EIR')

ggplot() +
  geom_sf(data = Senegambia_neighbors, fill = "cornsilk2", color = "cornsilk3") +
  geom_sf(data = Senegambia[Senegambia$ID_0 == 'GMB',], fill = "#72B000", color = "cornsilk3") +
  geom_sf(data = Senegambia[Senegambia$ID_0 == 'SEN',], fill = "#72B000", color = "cornsilk3") +
  geom_sf(data = plot_data, aes(fill = diff)) +
  theme_bw(base_size = 14) +
  # https://colorbrewer2.org/#type=diverging&scheme=RdBu&n=3
  scale_fill_gradient2(low = '#67a9cf', mid = "white", high = '#ef8a62', midpoint = 0) +
  scale_x_continuous(limits = c(-17.8, -11)) +
  scale_y_continuous(limits = c(12, 17)) +
  labs(x = '', y = '', fill = 'PfPR after \nweighted mixing') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "#daeff8", color = NA))

ggsave('./03_output/Senegambia_mixing_W.pdf', width = 6, height = 6)



