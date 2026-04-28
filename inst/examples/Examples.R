#### Example workflows for RankingDynamics
####
#### This script shows how to:
####   1. Convert raw score/abundance tables into ranked lists.
####   2. Run rank_analysis() on one ranked dataframe.
####   3. Run a permutation workflow to assess sensitivity to tied scores.
####   4. Summarize and plot results across datasets.
####
#### Expected input format:
####   - One row per entity, e.g. species, OTU, city, player.
####   - First column contains the entity identifier.
####   - Remaining columns contain scores/abundances at ordered time points.
####   - Larger values are treated as better/higher-ranked.
####   - Missing or non-positive values can be set to NA before ranking.

# -------------------------------------------------------------------------
# 0. Setup
# -------------------------------------------------------------------------

# If running this file inside the package source tree, source the package code.
# If the package is installed, you can replace these two lines with:
# library(RankingDynamics)
source("R/rank_analysis.R")
source("R/utils.R")

library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(1)

# -------------------------------------------------------------------------
# 1. Helper functions for examples
# -------------------------------------------------------------------------

# Generate an AR(1)-like score trajectory for one entity.
generate_ar1 <- function(n_time, rho = 1, sd = 1) {
  x <- numeric(n_time)
  x[1] <- rnorm(1, sd = sd)

  for (tt in 2:n_time) {
    x[tt] <- rho * x[tt - 1] + rnorm(1, sd = sd)
  }

  x
}

# Simulate a score table in the format expected by get_ranked_data().
simulate_drift_scores <- function(n_entities = 700, n_time = 100, rho = 1) {
  scores <- replicate(
    n = n_entities,
    expr = generate_ar1(n_time = n_time, rho = rho)
  )

  scores <- t(scores)
  scores[scores <= 0] <- NA

  colnames(scores) <- paste0("time_", seq_len(n_time))

  data.frame(
    entity = paste0("sp_", seq_len(n_entities)),
    scores,
    check.names = FALSE
  )
}

# Read a score table and convert non-positive scores to NA.
read_score_table <- function(path) {
  if (!file.exists(path)) {
    warning("File not found: ", path)
    return(NULL)
  }

  x <- read.csv(path, check.names = FALSE)
  x[-1][x[-1] <= 0] <- NA
  x
}

# Run the standard single-analysis workflow.
run_single_rank_analysis <- function(wide_data, N0 = NULL, verbose = FALSE) {
  ranked_data <- get_ranked_data(wide_data, N0 = N0)
  rank_analysis(ranked_data, verbose = verbose)
}

# Run permutation-based sensitivity analysis for tied scores.
#
# This does not change rank_analysis() itself. Instead, it creates many
# slightly jittered versions of the raw score table, converts each to ranks,
# and runs rank_analysis() on each ranked dataframe.
run_permuted_rank_analysis <- function(wide_data,
                                       n_perm = 100,
                                       N0 = NULL,
                                       epsilon = NULL,
                                       verbose = FALSE) {
  permuted_scores <- permute_dataframe(
    wide_data = wide_data,
    N = n_perm,
    epsilon = epsilon
  )

  ranked_perms <- get_ranked_data(permuted_scores, N0 = N0)

  lapply(
    ranked_perms,
    rank_analysis,
    verbose = verbose
  )
}

# Extract one-row summaries from rank_analysis() output.
summarise_rank_result <- function(result, dataset) {
  data.frame(
    dataset = dataset,
    p = result$p,
    mean_rank_flux = result$mean_rank_flux,
    mean_rank_turnover = result$mean_rank_turnover,
    tau = result$ptau_star,
    nu = result$pnu_star,
    tau_r = result$tau_r,
    nu_r = result$nu_r,
    Wlevy = result$Wlevy,
    Wdiff = result$Wdiff,
    Wrepl = result$Wrepl,
    fit_quality = result$fit_quality,
    check_1 = result$check_1,
    check_2 = result$check_2,
    check_3 = result$check_3,
    stringsAsFactors = FALSE
  )
}

# Summarize the distribution of model outputs across permutations.
summarise_permuted_results <- function(results, dataset) {
  summary_df <- bind_rows(
    lapply(results, summarise_rank_result, dataset = dataset),
    .id = "permutation"
  )

  summary_df %>%
    group_by(dataset) %>%
    summarise(
      n_perm = n(),
      across(
        c(p, mean_rank_flux, mean_rank_turnover, tau, nu, tau_r, nu_r,
          Wlevy, Wdiff, Wrepl, check_1, check_2, check_3),
        list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)),
        .names = "{.col}_{.fn}"
      ),
      n_failed = sum(fit_quality != "ok", na.rm = TRUE),
      .groups = "drop"
    )
}

# Create time-series plotting data from rank_analysis() output.
make_time_series_df <- function(result, dataset) {
  tibble(
    dataset = dataset,
    time_index = seq_len(result$Time),
    timefrac = result$timefrac,
    rank_turnover = result$rank_turnover,
    rank_flux = c(NA, result$rank_flux),
    rank_inertia = c(NA, result$rank_inertia)
  )
}

# -------------------------------------------------------------------------
# 2. Simulated drift example
# -------------------------------------------------------------------------

wide_data_drift <- simulate_drift_scores(
  n_entities = 700,
  n_time = 100,
  rho = 1
)

results_drift <- run_single_rank_analysis(wide_data_drift)

# Optional: permutation sensitivity analysis.
# For a simulated continuous process, ties are usually rare, so this may show
# little variation. It is more useful for ecological abundance tables with many
# tied counts.
results_drift_perm <- run_permuted_rank_analysis(
  wide_data = wide_data_drift,
  n_perm = 100,
  N0 = NULL
)

drift_perm_summary <- summarise_permuted_results(
  results = results_drift_perm,
  dataset = "Drift"
)

print(drift_perm_summary)

# -------------------------------------------------------------------------
# 3. Empirical examples
# -------------------------------------------------------------------------

# These paths assume your original project structure. Missing files are skipped.
data_paths <- c(
  BHC = "data/Pitcherplant/pitcher_otus_bhc.csv",
  BVBA = "data/Pitcherplant/pitcher_otus_bvba.csv",
  BCI = "data/BCI/bci_dbh100.csv"
)

wide_empirical <- lapply(path.expand(data_paths), read_score_table)
wide_empirical <- wide_empirical[!vapply(wide_empirical, is.null, logical(1))]

wide_empirical <- lapply(wide_empirical, function(df) {
  bad <- is.na(names(df)) | names(df) == ""
  names(df)[bad] <- paste0("unnamed_", which(bad))
  names(df) <- make.names(names(df), unique = TRUE)
  df
})

results_empirical <- lapply(
  wide_empirical,
  run_single_rank_analysis,
  N0 = NULL,
  verbose = FALSE
)

names(results_empirical) <- c("BHC_pitcher_microbes", "BVBA_pitcher_microbes", "BCI_trees")

all_results <- c(
  list(Drift = results_drift),
  results_empirical
)
# One-row summary per dataset.
result_summary <- bind_rows(
  Map(summarise_rank_result, all_results, names(all_results))
)

print(result_summary)

# -------------------------------------------------------------------------
# 4. Plot time-dependent quantities
# -------------------------------------------------------------------------

time_series_df <- bind_rows(
  Map(make_time_series_df, all_results, names(all_results))
)

plot_turnover <- ggplot(time_series_df, aes(x = timefrac, y = rank_turnover, color = dataset)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Rank turnover",
    x = "Normalized time",
    y = expression(o[t] == N[t] / N[0])
  ) +
  theme_minimal()

plot_flux <- ggplot(time_series_df, aes(x = timefrac, y = rank_flux, color = dataset)) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  labs(
    title = "Rank flux",
    x = "Normalized time",
    y = expression(F[t])
  ) +
  theme_minimal()

plot_inertia <- ggplot(time_series_df, aes(x = timefrac, y = rank_inertia, color = dataset)) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  labs(
    title = "Rank inertia",
    x = "Normalized lag",
    y = expression(S^{"++"})
  ) +
  theme_minimal()

# -------------------------------------------------------------------------
# 5. Plot model regimes in rescaled nu/tau space
# -------------------------------------------------------------------------

curve_df <- tibble(
  nu_r = 10^seq(-3, 3, length.out = 1000),
  tau_r = 1 / nu_r
)

plot_regime <- result_summary %>%
  filter(is.finite(nu_r), is.finite(tau_r), nu_r > 0, tau_r > 0) %>%
  ggplot(aes(x = nu_r, y = tau_r, color = dataset)) +
  geom_line(
    data = curve_df,
    aes(x = nu_r, y = tau_r),
    inherit.aes = FALSE,
    linetype = "dashed"
  ) +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Rescaled model parameters",
    x = expression(nu[r]),
    y = expression(tau[r])
  ) +
  theme_minimal()

combined_plot <- (plot_turnover + plot_flux) / (plot_inertia + plot_regime) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(combined_plot)

# Save if desired.
# ggsave("ranking_dynamics_examples.pdf", combined_plot, width = 10, height = 10)

# -------------------------------------------------------------------------
# 6. Example: permutation uncertainty for one empirical dataset
# -------------------------------------------------------------------------

# Run this block if you want a full tie-sensitivity workflow.
# This can be slow for large abundance tables.
# Single non-permuted run
bci_single <- run_single_rank_analysis(
  wide_data = wide_empirical[[3]],
  N0 = NULL,
  verbose = FALSE
)

# Summarise single result in the same format
bci_single_summary <- data.frame(
  dataset = "BCI",
  analysis = "single",
  Wlevy_mean = bci_single$Wlevy,
  Wdiff_mean = bci_single$Wdiff,
  Wrepl_mean = bci_single$Wrepl,
  tau_r_mean = bci_single$tau_r,
  nu_r_mean = bci_single$nu_r,
  mean_rank_flux_mean = bci_single$mean_rank_flux,
  mean_rank_turnover_mean = bci_single$mean_rank_turnover,
  stringsAsFactors = FALSE
)

# Permuted summary
bci_perms <- run_permuted_rank_analysis(
  wide_data = wide_empirical[[3]],
  n_perm = 500,
  N0 = NULL
)

bci_perm_summary <- summarise_permuted_results(
  results = bci_perms,
  dataset = "BCI"
)

bci_perm_summary$analysis <- "permuted"

# Compare
bci_compare <- dplyr::bind_rows(
  bci_single_summary,
  bci_perm_summary
)

# numbers should align quite closely IF random tiebreaks in a single run
# do not induce bias in the results.
print(bci_compare[,1:9])
