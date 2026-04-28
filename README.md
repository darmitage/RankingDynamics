
# RankingDynamics

An R package for quantifying temporal dynamics of ranked ecological, biological, and socioeconomic systems using rank flux, turnover, inertia, and model-derived parameters from displacement–replacement theory.

## Features

- Calculate rank flux (F) and turnover (o)
- Estimate diffusion, Lévy displacement, and replacement regimes
- Fit τ and ν parameters from empirical ranking trajectories
- Evaluate rank inertia and rank stability
- Perform permutation-based sensitivity analyses for tied rankings
- Reproduce major analyses from Iñiguez et al. (2022)

## Installation

```r
devtools::install_github("darmitage/RankingDynamics")

## Example

```r
# Load package
library(RankingDynamics)

# Convert abundance table into ranked trajectories
ranked_data <- get_ranked_data(wide_data)

# Run rank dynamics analysis
results <- rank_analysis(ranked_data)

# View fitted dynamical regime
results$Wlevy
results$Wdiff
results$Wrepl


## Citation

If you use this package, please cite:

> Iñiguez G, Pineda C, Gershenson C, Barabási A-L, Karsai M, Morales AJ.  
> **Dynamics of ranking.**  
> *Nature Communications* (2022) 13:1646.  
> https://doi.org/10.1038/s41467-022-29256-x