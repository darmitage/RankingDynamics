# RankingDynamics

An R package for quantifying temporal dynamics in ranked ecological, biological, and socioeconomic systems using rank flux, turnover, inertia, and displacement–replacement theory.

`RankingDynamics` provides tools to transform longitudinal abundance or score data into ranked trajectories, estimate empirical stability metrics, and fit the mechanistic framework described by Iñiguez et al. (2022) to identify whether system dynamics are dominated by diffusion, Lévy-like displacement, or replacement.

---

## Overview

Ranking systems are common across nature and society, from microbial communities and forest composition to economic systems and social hierarchies. This package allows users to:

- Convert temporal abundance or score tables into ranked trajectories
- Quantify:
  - Mean rank flux (`F`)
  - Rank turnover (`o`)
  - Rank inertia (`S++`)
  - Rank-change probabilities
- Fit model parameters:
  - `τ` (displacement)
  - `ν` (replacement)
- Estimate dynamical regimes:
  - `Wlevy`
  - `Wdiff`
  - `Wrepl`
- Assess uncertainty from tied ranks using permutation approaches
- Reproduce key analyses from published rank-dynamics theory

---

## Features

### Core functions

- `get_ranked_data()` — Convert wide-format score/abundance tables into rank trajectories
- `rank_analysis()` — Calculate rank metrics and fit displacement–replacement model parameters
- `permute_dataframe()` — Randomize tied observations for uncertainty estimation
- `run_permuted_rank_analysis()` — Perform repeated analyses across permutations
- `summarise_permuted_results()` — Summarize uncertainty distributions

### Supported applications

- Community ecology
- Species abundance dynamics
- Microbiome turnover
- Forest succession
- Evolutionary systems
- Socioeconomic rankings
- General ranked longitudinal datasets

---

## Installation

Install the development version from GitHub:

```r
devtools::install_github("darmitage/RankingDynamics")
```

Load the package:

```r
library(RankingDynamics)
```

---

## Basic Example

### Single dataset analysis

```r
# Load data
ranked_data <- get_ranked_data(wide_data)

# Analyze rank dynamics
results <- rank_analysis(ranked_data)

# Examine dynamical regime
results$Wlevy
results$Wdiff
results$Wrepl
```

---

## Permutation-Based Sensitivity Analysis

For datasets with many tied ranks:

```r
perm_results <- run_permuted_rank_analysis(
  wide_data = wide_data,
  n_perm = 500
)

summary <- summarise_permuted_results(
  results = perm_results,
  dataset = "ExampleDataset"
)

print(summary)
```

This workflow estimates how sensitive inferred ranking dynamics are to tie structure.

---

## Included Example Data

The package includes example datasets from:

- **BCI tropical forest dynamics**
- **Pitcher plant microbial communities**
- **Simulated drift communities**

Example workflows are available in:

```txt
inst/examples/Examples.R
```

---

## Repository Structure

```txt
RankingDynamics/
├── R/                  # Core package functions
├── data/               # Example datasets
│   ├── BCI/
│   └── Pitcherplant/
├── inst/examples/      # Reproducible example workflows
├── man/                # Documentation
└── README.md
```

---

## Scientific Basis

This package is based primarily on:

> Iñiguez G, Pineda C, Gershenson C, Barabási A-L, Karsai M, Morales AJ.  
> **Dynamics of ranking.**  
> *Nature Communications* (2022) 13:1646.  
> https://doi.org/10.1038/s41467-022-29256-x

### BibTeX

```bibtex
@article{iniguez2022ranking,
  title={Dynamics of ranking},
  author={Iñiguez, Gerardo and Pineda, Carlos and Gershenson, Carlos and Barabási, Albert-László and Karsai, Márton and Morales, Alfredo J.},
  journal={Nature Communications},
  volume={13},
  pages={1646},
  year={2022},
  doi={10.1038/s41467-022-29256-x}
}
```

### R citation

```r
citation("RankingDynamics")
```

---

## Future Directions

Planned improvements include:

- Native plotting functions for Figure S17-style visualizations
- Bayesian parameter fitting extensions
- Expanded null-model generation
- Greater compatibility with ecological time-series workflows
- CRAN submission

---

## Author

**David Armitage**  
Okinawa Institute of Science and Technology (OIST)

---

## License

MIT License
