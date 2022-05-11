# libraries
library(tidyverse)
library(patchwork)
library(ICDMM)

# define parameters [Marshall et al. 2018]
a <- 1.91
log_p <- 4.29
t <- 1.22

# destination population matrix
Nj <- matrix(data = c(1000, 1000, 1000, 1000,
                      500, 500, 500, 500,
                      500, 500, 500, 500,
                      300, 300, 300, 300),
             nrow = 4,
             ncol = 4)

# distance matrix
d <- matrix(data = c(0, 1, 1, 1,
                     1, 0, 1, 1,
                     1, 1, 0, 1,
                     1, 1, 1, 0),
           nrow = 4,
           ncol = 4)

# calculate travel
totals <- (Nj ^ t) * (1 + (d / exp(log_p)))^(-a)

# standardize so that each location sums to 1
Pij <- totals / rowSums(totals)

# rebuild ICDMM
# odin::odin_package(getwd())
# click 'install and rebuild'

# run ICDMM
age_vector <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5,
                7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
year <- 365

EIR_vector <- c(1, 5, 10, 20)
n <- length(EIR_vector)

# isolated example
mixing_I <- diag(n)

# perfectly mixed example
mixing_M <- matrix(rep((1 / n), n * n),
                 nrow = n,
                 ncol = n)

# distance / population weighted
mixing_W <- (Pij)

# define function to run model
mix_model <- function(EIR_vector, n, mixing){
  out <- run_model_metapop(model = "odin_model_metapop",
                           het_brackets = 5,
                           age = age_vector,
                           init_EIR = EIR_vector,
                           init_ft = 0.4,
                           mix_h = mixing,
                           mix_v = ,
                           time = 5 * year,
                           # num_int = num_int,
                           # ITN_on = 1,
                           # ITN_interval = 3 * year,
                           irs_cov = rep(0, n),
                           itn_cov= rep(0, n))

  model <- cbind(out$t, out$prev2to10) %>% as_tibble()
  colnames(model) <- c('t', EIR_vector)
  model <- model %>%
    pivot_longer(cols = 2:(n + 1), names_to = 'EIR', values_to = 'prev2to10') %>%
    mutate(EIR_f = factor(EIR, levels = EIR_vector, labels = EIR_vector))

  ggplot(data = model) +
    geom_line(aes(x = t/365, y = prev2to10, group = EIR_f, color = EIR_f)) +
    scale_y_continuous(limits = c(0,0.5)) +
    labs(x = 'time (years)',
         y = 'PfPR 2-10',
         color = 'EIR') +
    theme_classic()

}

# run three scenarios
A <- mix_model(EIR_vector, n, mixing_I) + labs(title = 'Isolated')
B <- mix_model(EIR_vector, n, mixing_M) + labs(title = 'Perfectly mixed')
C <- mix_model(EIR_vector, n, mixing_W) + labs(title = 'Weighted')

# plot
A + B + C + plot_layout(guides = "collect", nrow=1) + plot_annotation(tag_levels = 'A')

# save
ggsave()

# to add:
lapply(EIR_vector, function(x){paste(':',x)})

