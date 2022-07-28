
# Specify HPC options ----------------------------------------------------------
library(didehpc)
setwd("M:/Hillary/E8_spatial_mixing/")

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "htopazia")

# transfer the new malariasimulation folder manually to contexts or delete and re-install using conan
# remotes::install_github('mrc-ide/malariasimulation@dev', force=T)
src <- conan::conan_sources("github::mrc-ide/malariasimulation@bug/metapop_rowsums")

ctx <- context::context_save(path = "Q:/contexts",
                             sources = c('M:/Hillary/E8_spatial_mixing/function_metapop_sims.R'),
                             packages = c("dplyr", "malariasimulation"),
                             package_sources = src)

share <- didehpc::path_mapping("malaria", "M:", "//fi--didef3.dide.ic.ac.uk/malaria", "M:")
config <- didehpc::didehpc_config(shares = share,
                                  use_rrq = FALSE,
                                  cores = 1,
                                  cluster = "fi--didemrchnb", # fi--dideclusthn OR fi--didemrchnb
                                  parallel = FALSE)

obj <- didehpc::queue_didehpc(ctx, config = config)


# Set up your job --------------------------------------------------------------

# run weighted matrix
t <- obj$enqueue_bulk('W', metapop_sims)

# run isolated matrix
t <- obj$enqueue_bulk('I', metapop_sims)


# Read in results --------------------------------------------------------------
output1 <- readRDS('M:/Hillary/E8_spatial_mixing/metapop_run_I.rds') %>% mutate(model = 'isolated')
output2 <- readRDS('M:/Hillary/E8_spatial_mixing/metapop_run_W.rds') %>% mutate(model = 'mixed')

master <- readRDS('M:/Hillary/E8_spatial_mixing/E8_master.rds') %>% st_drop_geometry()

# combine results
output <- output1 %>% rbind(output2) %>%
  filter(year %in% c(7, 8, 9)) %>%
  group_by(country, admin1, model) %>%
  summarize(pfpr2_10 = mean(pfpr2_10)) %>%
  left_join(master %>% dplyr::select(COUNTRY, NAME_1, PfPR_rmean), by = c("country" = "COUNTRY", "admin1" = "NAME_1")) %>%
  ungroup() %>%
  mutate(PfPR_rmean = ifelse(is.na(PfPR_rmean), 0, PfPR_rmean),
         PRdiff = PfPR_rmean - pfpr2_10,
         index = row_number())


# plot difference between MAP PR vs. model output PR
ggplot(data = output, aes(x = index, y = PRdiff)) +
  geom_hline(aes(yintercept = 0), alpha = 0.3) +
  geom_point(aes(color = model), show.legend = T) +
  labs(x = '',
       y = 'difference in PfPR 2-10',
       title = 'Target PfPR 2-10 vs. calibrated model PfPR 2-10',
       color = '') +
  theme_classic() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

ggsave('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing/03_output/calibrate1.pdf', width = 5, height = 3)


# plot MAP PR vs. model output PR
ggplot(data = output, aes(x = PfPR_rmean, y = pfpr2_10)) +
  geom_abline(aes(intercept = 0, slope = 1), alpha = 0.3) +
  geom_point(aes(color = country), alpha = 0.7, show.legend = T) +
  facet_wrap(~model) +
  labs(x = 'MAP calibration PfPR 2-10',
       y = 'malariasimulation PfPR 2-10',
       color = '') +
  theme_classic()

ggsave('C:/Users/htopazia/OneDrive - Imperial College London/Github/E8_spatial_mixing/03_output/calibrate2.pdf', width = 5, height = 3)

