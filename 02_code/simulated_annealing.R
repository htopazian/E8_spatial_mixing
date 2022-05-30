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
