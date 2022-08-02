library(tidyverse)
library(readxl)
library(LaCroixColoR)

# print palette colors
lacroix_palettes
lacroix_palette("Pamplemousse", type = "discrete")
"#EA7580" "#F6A1A5" "#F8CD9C" "#1BB6AF" "#088BBE" "#172869"
 "1"       "3"       "2"       "5"       "4"       "6"

funding <- readxl::read_excel("./01_data/WMR2021_Annex5C.xlsx", sheet = "E8") %>%
  dplyr::select(Country, Year, donor_aid_reported, Gov, country_aid_reported) %>%
  pivot_longer(cols = donor_aid_reported:country_aid_reported,
               names_to = 'var', values_to = 'dollars')

A <- ggplot(data = funding %>% filter(var != 'country_aid_reported')) +
  geom_col(aes(x = Country, y = dollars / 1000000, fill = var), alpha = 0.8) +
  facet_wrap(~ Year, nrow = 3) +
  scale_fill_manual(values = c("tomato", "#088BBE"), labels = c('International funding', 'Government funding')) +
  labs(y = 'US dollars (millions, 2020)',
       x = '',
       fill = '',
       title = 'Donor-reported funding') +
  theme_classic()

B <- ggplot(data = funding %>% filter(var != 'country_aid_reported')) +
  geom_col(aes(x = Country, y = dollars / 1000000, fill = var), alpha = 0.8, show.legend = F) +
  facet_wrap(~ Year, nrow = 3) +
  scale_fill_manual(values = c("tomato", "#088BBE")) +
  labs(y = 'US dollars (millions, 2020)',
       x = '',
       fill = '',
       title = 'Country-reported funding') +
  theme_classic()


# consolidate plots
mapplot <- A + B + plot_layout(guides = "collect", nrow = 1) +
  plot_annotation(tag_levels = 'A')

ggsave('./03_output/funding.pdf', plot = mapplot, width = 14, height = 6)



# map flow
# E8 map -----------------------------------------------------------------------
countries <- readRDS('./03_output/countries.rds')
E8_sf <- countries %>% filter(ID_0 %in% E8) %>% mutate(COUNTRY = ifelse(COUNTRY == 'Swaziland', 'Eswatini', COUNTRY))
weighted_centroids <- readRDS('./03_output/E8_weighted_centroids.rds') %>% st_as_sf(crs = crs(countries))
master <- readRDS('./03_output/E8_master.rds')

# create point for legend
x <- c(1); y <- c(1); point <- tibble(x, y)

# plot flow between centroids
data <- readRDS('./03_output/E8_master.rds')
matrix <- readRDS('M:/Hillary/E8_spatial_mixing/E8_mixing_W.rds')

test <- as.data.frame(matrix)
test <- test %>% pivot_longer(cols = everything(), names_to = 'destination', values_to = 'flow')

test$country_origin <- unlist(lapply(master$ID_0, function(x){rep(x, nrow(master))}))
test$admin1_origin <- unlist(lapply(master$NAME_1, function(x){rep(x, nrow(master))}))
test$country_destination <- rep(master$ID_0, nrow(master))
test$admin1_destination <- rep(master$NAME_1, nrow(master))

test2 <- test %>%
  left_join(weighted_centroids, by = c('country_origin' = 'ID_0', 'admin1_origin' = 'NAME_1')) %>%
  rename(centroid_origin = w_cent_geometry) %>%
  left_join(weighted_centroids, by = c('country_destination' = 'ID_0', 'admin1_destination' = 'NAME_1')) %>%
  rename(centroid_destination = w_cent_geometry) %>%
  mutate(long_origin = unlist(map(centroid_origin, 1)),
         lat_origin = unlist(map(centroid_origin, 2)),
         long_destination = unlist(map(centroid_destination, 1)),
         lat_destination = unlist(map(centroid_destination, 2)),
         d_lat = (lat_destination - lat_origin),
         d_long = (long_destination - long_origin)) %>%
  mutate(flow_f = case_when(flow == 0 ~ 'none',
                            flow > 0 & flow <= 0.001 ~ 'less movement',
                            flow > 0.001 & flow <= 0.1 ~ 'more movement',
                            flow > 0.1 ~ 'max movement'),
         flow_f = factor(flow_f)) %>%
  filter(flow > 0 & flow < 0.9) %>%
  left_join(E8_sf %>% dplyr::select(ID_0, COUNTRY), by = c('country_origin' = 'ID_0'))

summary(test2$flow)
table(test2$flow_f)

library(patchwork)
library(grid) # needed for arrows argument

'A')

ggplot() +
  geom_sf(data = E8_sf, fill = "cornsilk2", color = "cornsilk3", size = 0.5) +
  geom_segment(data = test2 %>% filter(country_origin == country_destination),
               aes(x = long_origin, y = lat_origin, xend = long_destination, yend = lat_destination, color = flow),
               size = 1,
               alpha = 0.2,
               arrow = arrow(length = unit(0.05, 'cm'))) +
  # scale_alpha_continuous(range = c(.1, .7)) +
  scale_color_distiller(palette = 'Spectral') +
  theme_classic(base_size = 14) +
  facet_wrap(~ COUNTRY, scales = 'free') +
  # coord_sf(xlim = c(10, 51.5), ylim = c(-35, 6)) +
  labs(x = '', y = '', color = '', fill = '', color = '') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = NA, color = NA))

# consolidate plots
mapplot <- A + B + plot_layout(guides = "collect", nrow = 1) +
  plot_annotation(tag_levels =

summary(test2$flow)

beepr::beep(1)
ggsave('./03_output/flow.pdf', plot = C, width = 6, height = 6)



# internal mixing

plot_admin <- function(x) {

ggplot() +
  geom_sf(data = E8_sf %>% filter(COUNTRY == x), fill = "cornsilk2", color = "cornsilk3", size = 0.5) +
  geom_segment(data = test2 %>%
                 filter(country_origin == country_destination) %>%
                 filter(COUNTRY == x),
               aes(x = long_origin, y = lat_origin, xend = long_destination, yend = lat_destination, color = flow),
               size = .8,
               arrow = arrow(length = unit(0.1, 'cm')),
               alpha = 0.01) +
  scale_color_distiller(palette = 'Spectral', breaks = seq(0, 0.06, 0.01), limits = c(0, 0.06)) +
  theme_bw(base_size = 14) +
  labs(x = '', y = '', color = 'mixing', fill = '', color = '', title = x) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = NA, color = NA),
        panel.border = element_blank())

}

A <- plot_admin('Angola')
B <- plot_admin('Botswana')
C <- plot_admin('Eswatini')
D <- plot_admin('Mozambique')
E <- plot_admin('Namibia')
F <- plot_admin('South Africa')
G <- plot_admin('Zambia')
H <- plot_admin('Zimbabwe')


# consolidate plots
mapplot <- A + B + C + D + E + F + G + H + plot_layout(guides = "collect", nrow = 2, ncol = 4)

ggsave('./03_output/flow_admin1.pdf', plot = mapplot, width = 12, height = 8)


