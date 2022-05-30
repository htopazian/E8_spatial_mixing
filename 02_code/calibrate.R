# calibrate pfprs to EIRs in odin
# https://github.com/mrc-ide/cali/blob/main/R/calibrate.R
calibrate <- function(target, target_tt, summary_function, tolerance, interval = c(0.01, 2000) / 365, ...){
  stats::uniroot(objective,
                 target = target,
                 target_tt = target_tt,
                 summary_function = summary_function,
                 tolerance = tolerance,
                 interval = interval,
                 ...)
}

objective <- function(x, target, target_tt, summary_function, tolerance){
  message("Trying EIR: ", signif(x, 3))

  raw_output <- ICDMM::run_model(model = "odin_model",
                                 init_EIR = x,
                                 init_ft = 0.4,
                                 time = target_tt)

  raw_output <- as_tibble(raw_output)

  target_variable <- raw_output$prev2to10
  difference <- (target_variable[target_tt] - target)
  message("Current difference: ", paste(signif(difference, 3), collapse = " "))
  # Adjust for specified tolerance
  difference[abs(difference) < tolerance] <- 0
  difference <- sum(difference)
  return(difference)
}

summary_pfpr_2_10 <- function(x){
  prev_2_10 <- x$prev2to10
  return(prev_2_10)
}


PR_EIR <- function(target, target_tt){
  set.seed(123)
  out <- calibrate(target = target,
                   target_tt = target_tt,
                   summary_function = summary_pfpr_2_10,
                   tolerance = 0.02,
                   interval = c(.0001, 300))

  root <- as_tibble(out) %>% dplyr::select(root)

  root

}
