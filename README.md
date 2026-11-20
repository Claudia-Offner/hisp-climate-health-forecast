# Associations between extreme weather and healthcare use in the Indian Himalayas

This repository applies a Distributed Lag Non-linear Model (DLNM) time series approach to characterise the effects of weather on hospital attendance in Himachal Pradesh, India. Code is provided for transparency and future reference. Some of the data are not publicly available; please direct data enquiries to the authors.

## Code directory

The full analysis is run from `main.R`. Shared dependencies are in `utilities.R`.

- `1_parameter_selection.R`: assesses and selects parameters for the DLNM models
- `2_statistical_analysis.R`: runs the main analysis, including meta-models, predictor tests and subgroup tests
- `3_figures_and_tables.R`: compiles the main and supplementary figures and tables

## Outputs

- `outputs/`: model outputs from `1_parameter_selection.R` and `2_statistical_analysis.R`
- `results/`: main and supplementary figures and tables from `3_figures_and_tables.R`, plus sensitivity analyses from `1_parameter_selection.R`