#### Main code for calculating Rank Flux model paramters (see below)

#' @title Analyze Rank Dynamics in a Time Series
#'
#' @description This function performs an analysis of rank dynamics from a time-series
#' dataset. The dataset should include ranks of entities across multiple time points.
#' The function computes metrics like rank flux, probability of leaving ranks,
#' probability of rank changes, rank turnover, and rank inertia. It also solves
#' transcendental equations for rank flux and openness metrics, performing checks for validity.
#'
#' @param rank_data A data.frame where the first column contains entities (e.g., species)
#' and subsequent columns represent rankings or scores at different time points.
#'
#' @return A list containing the computed metrics and results:
#' - `rank_data`: Processed rank data.
#' - `Wlevy`, `Wdiff`, `Wrepl`: Metrics for ranking dynamics.
#' - `pnu_star`, `ptau_star`: Solutions to transcendental equations.
#' - `tau_r`, `nu_r`: Time-related metrics.
#' - `p`: Probability of remaining ranked.
#' - `Spp`: Species persistence.
#' - `rank_inertia`: Rank inertia over time.
#' - `Nt`: Turnover-adjusted rank counts.
#' - `Time`: Number of time points.
#' - `rank_turnover`: Rank turnover at each time point.
#' - `mean_rank_turnover`: Average rank turnover.
#' - `rank_flux`: Rank flux for transitions between time points.
#' - `mean_rank_flux`: Average rank flux.
#' - `leave_probabilities`: Probabilities of leaving ranks.
#' - `probabilities_rank_change`: Probabilities of rank changes.
#' - `timefrac`: Normalized time fraction for plotting.
#' - `check_1`, `check_2`, `check_3`: Results of validity checks for calculations.
#'
#' @importFrom stats optimize
#' @import dplyr

#' @export

rank_analysis <- function(rank_data) {

  rank_data <- rank_data[,-1] # remove first column
  N0 <- length(na.omit(rank_data[,1])) # estimate N0 from data

  # Rank Flux Calculation (values range between zero and one, indicating a )
  # measure of openness or closedness
  # NEED TO MODIFY FOR CLOSED LISTS
  rank_flux_list <- list()
  for (i in 2:ncol(rank_data)) {
    top_t_minus_1 <- which(rank_data[, i-1] <= N0)
    top_t <- which(rank_data[, i] <= N0)
    leaves <- setdiff(top_t_minus_1, top_t)
    flux <- length(leaves) / N0
    rank_flux_list[[paste0(colnames(rank_data)[i-1], " to ", colnames(rank_data)[i])]] <- flux
  }
  rank_flux <- unlist(rank_flux_list)
  mean_rank_flux <- mean(rank_flux)

  # Probability of Leaving Ranking
  leave_counts <- rep(0, N0)
  total_counts <- rep(0, N0)

  for (i in 2:ncol(rank_data)) {
    rank_t_minus_1 <- rank_data[, i - 1]
    rank_t <- rank_data[, i]
    left_ranking <- is.na(rank_t)
    for (R in 1:N0) {
      leave_counts[R] <- leave_counts[R] + sum(left_ranking & rank_t_minus_1 == R, na.rm = TRUE)
      total_counts[R] <- total_counts[R] + sum(rank_t_minus_1 == R, na.rm = TRUE)
    }
  }
  leave_probabilities <- leave_counts / total_counts

  # Probability of Rank Change
  change_counts <- rep(0, N0)

  for (i in 1:(ncol(rank_data) - 1)) {
    rank_t <- rank_data[, i]
    rank_t_plus_1 <- rank_data[, i + 1]
    for (R in 1:N0) {
      entity_t <- rownames(rank_data)[which(rank_t == R)]
      entity_t_plus_1 <- rownames(rank_data)[which(rank_t_plus_1 == R)]
      if (length(entity_t) == 1 && length(entity_t_plus_1) == 1) {
        if (entity_t != entity_t_plus_1) {
          change_counts[R] <- change_counts[R] + 1
        }
        total_counts[R] <- total_counts[R] + 1
      }
    }
  }
  probabilities_rank_change <- change_counts / total_counts

  # Rank Turnover
  unique_elements_seen <- list()
  rank_turnover <- numeric(ncol(rank_data))

  for (t in 1:ncol(rank_data)) {
    elements_in_t <- rownames(rank_data)[!is.na(rank_data[, t])]
    if (t == 1) {
      unique_elements_seen[[t]] <- elements_in_t
    } else {
      unique_elements_seen[[t]] <- unique(c(unique_elements_seen[[t - 1]], elements_in_t))
    }
    rank_turnover[t] <- length(unique_elements_seen[[t]]) / N0
  }

  # Rank Inertia
  top_half_threshold <- 0.5 * N0
  rank_inertia <- numeric(ncol(rank_data) - 1)
  initial_top_half <- rownames(rank_data)[rank_data[, 1] <= top_half_threshold & !is.na(rank_data[, 1])]

  for (t in 2:ncol(rank_data)) {
    top_half_in_t <- rownames(rank_data)[rank_data[, t] <= top_half_threshold & !is.na(rank_data[, t])]
    remaining_in_top_half <- sum(initial_top_half %in% top_half_in_t)
    rank_inertia[t - 1] <- remaining_in_top_half / length(initial_top_half)
  }

  mean_rank_turnover <- (tail(rank_turnover, 1) - rank_turnover[1]) / length(rank_turnover)
  Nt <- rank_turnover * N0
  Time <- ncol(rank_data)
  N <- tail(Nt, 1)
  p <- N0 / N

  # Transcendental Equation Calculations (using base R functions)
  flux <- mean_rank_flux
  open_deriv <- mean_rank_turnover
  p0 <- p
  Splusplus <- rank_inertia[1]

  ptau_func <- function(pnu) {
    pnu * (pnu - open_deriv) / (p0 * open_deriv - pnu)
  }
  exp_func <- function(pnu) {
    exp(-ptau_func(pnu))
  }
  pnu_func <- function(pnu) {
    log((p0 + (1 - p0) * exp_func(pnu)) / (1 - flux))
  }

  pnu_lambda <- function(pnu) {
    abs(pnu - pnu_func(pnu))
  }

  if(flux == 0 & open_deriv == 0){

    pnu_star <- 0

    ptau_star <- log((1-0.5*p0)/(Splusplus-0.5*p0))

  } else {

    pnu_res <- optimize(pnu_lambda,
                        interval = c(p0 * open_deriv, open_deriv),
                        tol = .Machine$double.eps^0.9)

    pnu_star <- pnu_res$minimum

    ptau_star <- ifelse(ptau_func(pnu_star) > 1, 1, ptau_func(pnu_star))
  }
  # Final Calculations
  tau_r <- ptau_star / (p0 * (1 - p0) * open_deriv)
  nu_r <- (pnu_star - p0 * open_deriv) / open_deriv
  Wlevy <- exp(-pnu_star) * (1 - exp(-ptau_star))
  Wdiff <- exp(-pnu_star) * exp(-ptau_star)
  Wrepl <- 1 - exp(-pnu_star)

  # Condition checks and warnings
  check_1 <- abs(round(pnu_star * ((pnu_star + ptau_star) / (pnu_star + p0 * ptau_star)), 12) - open_deriv)
  check_2 <- abs(round(1 - exp(-pnu_star) * (p0 + (1 - p0) * exp(-ptau_star)), 12) - flux)
  check_3 <- Wlevy + Wdiff + Wrepl

  # Print the results of the checks
  cat("Check 1: open_deriv - f(ptau_star, pnu_star) =", check_1, "<-should be near zero \n")
  cat("Check 2: flux_difference - f(ptau_star, pnu_star)  =", check_2, "<-should be near zero \n")
  cat("Check 3: (Wlevy + Wdiff + Wrepl difference) =", check_3, "<-should be 1 \n")

  timefrac = (1:Time) / Time
  # Return results
  list(
    rank_data = rank_data,
    Wlevy = Wlevy,
    Wdiff = Wdiff,
    Wrepl = Wrepl,
    pnu_star = pnu_star,
    ptau_star = ptau_star,
    tau_r = tau_r,
    nu_r = nu_r,
    p = p,
    Spp = Splusplus,
    rank_inertia = rank_inertia,
    Nt = Nt,
    Time = Time,
    rank_turnover = rank_turnover,
    mean_rank_turnover = mean_rank_turnover,
    rank_flux = rank_flux,
    mean_rank_flux = mean_rank_flux,
    rank_inertia = rank_inertia,
    leave_probabilities = leave_probabilities,
    probabilities_rank_change = probabilities_rank_change,
    timefrac = timefrac
  )
}
