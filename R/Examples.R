
#### This script loads and calculates ranking dynamics paramters for
#### four different datasets:
# 1. A simulated community undergoing only drift
# 2. Microbial abundances from a yearlong, uneven time series of
#     Darlingtonia californica pitcher plant goo
# 3. Microbial abundances from an yearlong, uneven series of
#     Darlingtonia californica pitcher plant goo from a different site
# 4. 50-ha plot data from BCI (all stems larger than 10 cm dbh)

# Load required libraries
source("R/rank_analysis.R")
source("R/utils.R")
library(dplyr)
library(ggplot2)
library(patchwork)

# This function generates a community of N species whose populations are subject
# to pure drift over T generations. We can set the rho value to change
# autocorrelative behavior of random walk

# Function to generate an AR(1) time series
generate_ar1 <- function(T, rho) {
  # Initialize the time series
  series <- numeric(T)
  # Initial value
  series[1] <- rnorm(1)
  # Generate subsequent values with AR(1) process
  for (t in 2:T) {
    series[t] <- rho * series[t-1] + rnorm(1)
  }
  return(series)
}

# Set parameters
N <- 700        # Number of species
T <- 100       # Number of time points
rho <- 1     # Autocorrelation coefficient (close to 1 means high autocorrelation)


population_abundances <- matrix(nrow = N, ncol = T) # initialise matrix

# simulate dynamics
for (i in 1:N) {
  population_abundances[i, ] <- generate_ar1(T, rho)
}
population_abundances[1:5, 1:10]  # First 5 species and first 10 time points
population_abundances[population_abundances <= 0] <- NA

wide_data <- data.frame(species = 1:N, population_abundances)

ranked_data_drift <- get_ranked_data(wide_data)
results_drift_ac <- rank_analysis(ranked_data_drift)

### Do the same thing with the pitcher plant microbial OTU abundance tables
### from blackhawk creek site (Plumas CO, CA)
wide_data <- read.csv("data/pitcher_otus_bhc.csv")
wide_data[wide_data <= 0] <- NA
ranked_data_bhc <- get_ranked_data(wide_data)
results_bhc <- rank_analysis(ranked_data_bhc)

### Do the same thing with the pitcher plant microbial OTU abundance tables
### from BVBA site (Plumas CO, CA)
wide_data <- read.csv("data/pitcher_otus_bvba.csv")
wide_data[wide_data <= 0] <- NA
ranked_data_bvba <- get_ranked_data(wide_data)
results_bvba <- rank_analysis(ranked_data_bvba)

### Do the same thing with the pitcher plant OTU relative abundance tables
### from 50-ha plot tree data at Barro Colorado Island (DBH > 100 mm)
wide_data <- read.csv("data/bci_dbh100.csv")
wide_data[wide_data <= 0] <- NA
ranked_data_bci <- get_ranked_data(wide_data)
results_bci <- rank_analysis(ranked_data_bci)

#### Use GGPLOT2 to plot the results!

plot1 <- ggplot() +
  geom_line(aes(x = results_drift_ac$timefrac, y = results_drift_ac$rank_turnover, color = "Drift AC"), size = 1.2) +
  geom_line(aes(x = results_bhc$timefrac, y = results_bhc$rank_turnover, color = "BHC"), size = 1.2) +
  geom_line(aes(x = results_bvba$timefrac, y = results_bvba$rank_turnover, color = "BVBA"), size = 1.2) +
  geom_line(aes(x = results_bci$timefrac, y = results_bci$rank_turnover, color = "BCI"), size = 1.2) +
  labs(title = "Rank Turnover", x = "Time Fraction", y = "Rank Turnover") +
  theme_minimal()

plot2 <- ggplot() +
  geom_line(aes(x = results_drift_ac$timefrac[-1], y = results_drift_ac$rank_flux, color = "Drift AC"), size = 1.2) +
  geom_line(aes(x = results_bhc$timefrac[-1], y = results_bhc$rank_flux, color = "BHC"), size = 1.2) +
  geom_line(aes(x = results_bvba$timefrac[-1], y = results_bvba$rank_flux, color = "BVBA"), size = 1.2) +
  geom_line(aes(x = results_bci$timefrac[-1], y = results_bci$rank_flux, color = "BCI"), size = 1.2) +
  labs(title = "Rank Flux", x = "Time Fraction", y = "Rank Flux") +
  theme_minimal()

plot3 <- ggplot() +
  geom_line(aes(x = results_drift_ac$timefrac[-1], y = results_drift_ac$rank_inertia, color = "Drift AC"), size = 1.2) +
  geom_line(aes(x = results_bhc$timefrac[-1], y = results_bhc$rank_inertia, color = "BHC"), size = 1.2) +
  geom_line(aes(x = results_bvba$timefrac[-1], y = results_bvba$rank_inertia, color = "BVBA"), size = 1.2) +
  geom_line(aes(x = results_bci$timefrac[-1], y = results_bci$rank_inertia, color = "BCI"), size = 1.2) +
  labs(title = "Rank Inertia", x = "Time Fraction", y = "Rank Inertia") +
  theme_minimal()

plot4 <- ggplot() +
  geom_line(aes(x = log(x), y = log(y)), linetype = "dashed") +
  geom_point(aes(x = log(results_drift_ac$nu_r), y = log(results_drift_ac$tau_r), color = "Drift AC"), size = 3) +
  geom_point(aes(x = log(results_bhc$nu_r), y = log(results_bhc$tau_r), color = "BHC"), size = 3) +
  geom_point(aes(x = log(results_bvba$nu_r), y = log(results_bvba$tau_r), color = "BVBA"), size = 3) +
  geom_point(aes(x = log(results_bci$nu_r), y = log(results_bci$tau_r), color = "BCI"), size = 3) +
  labs(title = "Nu ~ Tau Space", x = "log(nu)", y = "log(tau)") +
  theme_minimal()

# Combine plots with shared legend
combined_plot <- (plot1 + plot2) / (plot3 + plot4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Print the combined plot
print(combined_plot)

# ggsave("initial_ranking_results.pdf", plot = combined_plot, width = 10, height = 10)

