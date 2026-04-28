#' Analyze rank dynamics in a time series
#'
#' @description
#' Calculates empirical rank-dynamics summaries and fits the displacement /
#' replacement model described by Iniguez et al. The input should already be a
#' ranked data set: the first column contains entity identifiers and all remaining
#' columns contain ranks at successive time points. Ranks outside the retained
#' list should be coded as `NA`.
#'
#' @param rank_data A data frame. Column 1 contains entity identifiers; columns
#'   2 onward contain ranks through time.
#' @param verbose Logical. If `TRUE`, print model-fit diagnostics.
#' @param check_tol Numeric tolerance used to flag poor model fits.
#' @param inertia_threshold Numeric in `(0, 1)`. Fraction of `N0` defining the
#'   top region for rank inertia. The default, `0.5`, follows the supplement.
#' @param root_tol Numeric tolerance passed to `uniroot()`.
#'
#' @return A list containing empirical rank summaries, fitted model parameters,
#'   dynamical regime probabilities, and fit diagnostics.
#'
#' @export
rank_analysis <- function(rank_data,
                          verbose = TRUE,
                          check_tol = 1e-3,
                          inertia_threshold = 0.5,
                          root_tol = 1e-12) {
  validate_rank_data(rank_data, inertia_threshold)

  entity_names <- as.character(rank_data[[1]])
  ranks <- as.data.frame(rank_data[, -1, drop = FALSE])
  ranks[] <- lapply(ranks, as.numeric)

  n_time <- ncol(ranks)
  if (n_time < 2) {
    stop("rank_data must contain at least two rank/time columns after the entity column.", call. = FALSE)
  }

  # Ranking-list size. This assumes rank_data has already been trimmed to a
  # constant list size; get_ranked_data() does this by retaining the top N0 ranks.
  N0 <- sum(!is.na(ranks[[1]]))
  if (N0 <= 0) stop("The first rank column contains no non-NA ranks.", call. = FALSE)

  rank_flux <- calculate_rank_flux(ranks, N0)
  mean_rank_flux <- mean(rank_flux, na.rm = TRUE)

  leave_probabilities <- calculate_leave_probabilities(ranks, N0)
  probabilities_rank_change <- calculate_rank_change_probabilities(ranks, entity_names, N0)

  turnover <- calculate_rank_turnover(ranks, entity_names, N0)
  Nt <- turnover$Nt
  rank_turnover <- turnover$rank_turnover

  inertia <- calculate_rank_inertia(
    ranks = ranks,
    entity_names = entity_names,
    N0 = N0,
    c = inertia_threshold
  )

  p0 <- N0 / max(Nt)
  flux <- mean_rank_flux
  open_deriv <- turnover$mean_rank_turnover

  model <- fit_rank_model(
    p = p0,
    odot = open_deriv,
    F = flux,
    Splusplus = inertia$Splusplus,
    root_tol = root_tol
  )

  diagnostics <- calculate_model_diagnostics(
    p = p0,
    odot = open_deriv,
    F = flux,
    nu = model$pnu_star,
    tau = model$ptau_star,
    check_tol = check_tol
  )

  if (verbose) {
    message("Check 1 (turnover-rate match): ", signif(diagnostics$check_1, 5))
    message("Check 2 (flux match): ", signif(diagnostics$check_2, 5))
    message("Check 3 (Wlevy + Wdiff + Wrepl): ", signif(diagnostics$check_3, 5))
    message("Fit quality: ", diagnostics$fit_quality)
  }

  list(
    rank_data = ranks,
    entity_names = entity_names,
    N0 = N0,
    Time = n_time,
    p = p0,
    Nt = Nt,
    deltaNt = diff(Nt),
    rank_turnover = rank_turnover,
    mean_rank_turnover = open_deriv,
    rank_flux = rank_flux,
    mean_rank_flux = mean_rank_flux,
    leave_probabilities = leave_probabilities,
    probabilities_rank_change = probabilities_rank_change,
    rank_inertia = inertia$rank_inertia,
    rank_inertia_initial = inertia$rank_inertia_initial,
    Spp = inertia$Splusplus,
    pnu_star = model$pnu_star,
    ptau_star = model$ptau_star,
    nu_r = model$nu_r,
    tau_r = model$tau_r,
    Wlevy = diagnostics$Wlevy,
    Wdiff = diagnostics$Wdiff,
    Wrepl = diagnostics$Wrepl,
    check_1 = diagnostics$check_1,
    check_2 = diagnostics$check_2,
    check_3 = diagnostics$check_3,
    fit_quality = diagnostics$fit_quality,
    timefrac = seq_len(n_time) / n_time,
    transition_timefrac = seq(2, n_time) / n_time,
    inertia_timefrac = seq_len(n_time - 1) / n_time
  )
}

validate_rank_data <- function(rank_data, inertia_threshold) {
  if (!is.data.frame(rank_data)) {
    stop("rank_data must be a data frame.", call. = FALSE)
  }
  if (ncol(rank_data) < 3) {
    stop("rank_data must have one entity column and at least two rank/time columns.", call. = FALSE)
  }
  if (anyDuplicated(rank_data[[1]]) > 0) {
    stop("Entity identifiers in the first column must be unique.", call. = FALSE)
  }
  if (!is.numeric(inertia_threshold) || length(inertia_threshold) != 1 ||
      inertia_threshold <= 0 || inertia_threshold >= 1) {
    stop("inertia_threshold must be a single number between 0 and 1.", call. = FALSE)
  }
}

calculate_rank_flux <- function(ranks, N0) {
  vapply(2:ncol(ranks), function(t) {
    in_previous <- which(!is.na(ranks[[t - 1]]) & ranks[[t - 1]] <= N0)
    in_current <- which(!is.na(ranks[[t]]) & ranks[[t]] <= N0)
    length(setdiff(in_previous, in_current)) / N0
  }, numeric(1))
}

calculate_leave_probabilities <- function(ranks, N0) {
  leave_counts <- numeric(N0)
  total_counts <- numeric(N0)

  for (t in 2:ncol(ranks)) {
    previous <- ranks[[t - 1]]
    current <- ranks[[t]]

    for (rank in seq_len(N0)) {
      at_rank <- !is.na(previous) & previous == rank
      leave_counts[rank] <- leave_counts[rank] + sum(at_rank & is.na(current))
      total_counts[rank] <- total_counts[rank] + sum(at_rank)
    }
  }

  total_counts[total_counts == 0] <- NA_real_
  leave_counts / total_counts
}

calculate_rank_change_probabilities <- function(ranks, entity_names, N0) {
  change_counts <- numeric(N0)
  total_counts <- numeric(N0)

  for (t in 2:ncol(ranks)) {
    previous <- ranks[[t - 1]]
    current <- ranks[[t]]

    for (rank in seq_len(N0)) {
      previous_entity <- entity_names[!is.na(previous) & previous == rank]
      current_entity <- entity_names[!is.na(current) & current == rank]

      if (length(previous_entity) == 1 && length(current_entity) == 1) {
        change_counts[rank] <- change_counts[rank] + as.integer(previous_entity != current_entity)
        total_counts[rank] <- total_counts[rank] + 1
      }
    }
  }

  total_counts[total_counts == 0] <- NA_real_
  change_counts / total_counts
}

calculate_rank_turnover <- function(ranks, entity_names, N0) {
  seen <- character(0)
  Nt <- numeric(ncol(ranks))

  for (t in seq_len(ncol(ranks))) {
    current <- entity_names[!is.na(ranks[[t]]) & ranks[[t]] <= N0]
    seen <- union(seen, current)
    Nt[t] <- length(seen)
  }

  rank_turnover <- Nt / N0

  # Paper definition: odot = (o_{T-1} - o_0) / (T - 1), where there are
  # T observations and T - 1 transitions.
  mean_rank_turnover <-
    (tail(rank_turnover, 1) - rank_turnover[1]) / (length(rank_turnover) - 1)

  list(
    Nt = Nt,
    rank_turnover = rank_turnover,
    mean_rank_turnover = mean_rank_turnover
  )
}

calculate_rank_inertia <- function(ranks, entity_names, N0, c = 0.5) {
  top_cutoff <- c * N0
  n_time <- ncol(ranks)

  top_entities <- lapply(seq_len(n_time), function(t) {
    entity_names[!is.na(ranks[[t]]) & ranks[[t]] <= top_cutoff]
  })

  # Paper-style S++_lag: probability that an element in the top at t0 remains
  # in the top after 'lag' observations, averaged over all compatible t0.
  rank_inertia <- vapply(seq_len(n_time - 1), function(lag) {
    vals <- vapply(seq_len(n_time - lag), function(t0) {
      start_top <- top_entities[[t0]]
      end_top <- top_entities[[t0 + lag]]

      if (length(start_top) == 0) return(NA_real_)
      sum(start_top %in% end_top) / length(start_top)
    }, numeric(1))

    mean(vals, na.rm = TRUE)
  }, numeric(1))

  # Kept for comparison with the earlier implementation: retention relative only
  # to the first observed ranking.
  initial_top <- top_entities[[1]]
  rank_inertia_initial <- vapply(2:n_time, function(t) {
    if (length(initial_top) == 0) return(NA_real_)
    sum(initial_top %in% top_entities[[t]]) / length(initial_top)
  }, numeric(1))

  list(
    rank_inertia = rank_inertia,
    rank_inertia_initial = rank_inertia_initial,
    Splusplus = rank_inertia[1]
  )
}

fit_rank_model <- function(p, odot, F, Splusplus, root_tol = 1e-12) {
  if (isTRUE(all.equal(F, 0)) && isTRUE(all.equal(odot, 0))) {
    nu <- 0
    tau <- log((1 - 0.5 * p) / (Splusplus - 0.5 * p))

    return(list(
      pnu_star = nu,
      ptau_star = tau,
      nu_r = NA_real_,
      tau_r = NA_real_
    ))
  }

  root <- find_rank_model_root(p = p, odot = odot, F = F, tol = root_tol)
  nu <- root$nu_star
  tau <- root$tau_star

  list(
    pnu_star = nu,
    ptau_star = tau,
    nu_r = (nu - p * odot) / odot,
    tau_r = tau / (p * (1 - p) * odot)
  )
}

find_rank_model_root <- function(p, odot, F, tol = 1e-12, n_grid = 5000) {
  tau_of_nu <- function(nu) nu * (nu - odot) / (p * odot - nu)

  rhs_of_nu <- function(nu) {
    tau <- tau_of_nu(nu)
    log((p + (1 - p) * exp(-tau)) / (1 - F))
  }

  g <- function(nu) nu - rhs_of_nu(nu)

  pole <- p * odot
  hi <- odot
  eps <- .Machine$double.eps * 10

  grid <- seq(pole + eps, hi - eps, length.out = n_grid)
  g_values <- g(grid)
  keep <- is.finite(g_values)
  grid <- grid[keep]
  g_values <- g_values[keep]

  sign_change <- which(diff(sign(g_values)) != 0)
  if (length(sign_change) == 0) {
    stop("No model root found in the admissible interval.", call. = FALSE)
  }

  roots <- vapply(sign_change, function(i) {
    uniroot(g, interval = grid[i + 0:1], tol = tol)$root
  }, numeric(1))

  tau_roots <- tau_of_nu(roots)
  physical <- which(is.finite(tau_roots) & tau_roots > 0)

  if (length(physical) == 0) {
    stop("No physically admissible model root found.", call. = FALSE)
  }

  # Usually there is one physical root. If there are multiple, use the first.
  list(
    nu_star = roots[physical][1],
    tau_star = tau_roots[physical][1],
    extra_roots = roots[-physical]
  )
}

calculate_model_diagnostics <- function(p, odot, F, nu, tau, check_tol) {
  Wlevy <- exp(-nu) * (1 - exp(-tau))
  Wdiff <- exp(-nu - tau)
  Wrepl <- 1 - exp(-nu)

  if (is.na(odot) || isTRUE(all.equal(odot, 0))) {
    check_1 <- NA_real_
  } else {
    check_1 <- abs(nu * ((nu + tau) / (nu + p * tau)) - odot)
  }

  check_2 <- abs(1 - exp(-nu) * (p + (1 - p) * exp(-tau)) - F)
  check_3 <- Wlevy + Wdiff + Wrepl

  max_check <- max(check_1, check_2, na.rm = TRUE)
  fit_quality <- ifelse(max_check > check_tol, "fail", "ok")

  list(
    Wlevy = Wlevy,
    Wdiff = Wdiff,
    Wrepl = Wrepl,
    check_1 = check_1,
    check_2 = check_2,
    check_3 = check_3,
    fit_quality = fit_quality
  )
}
