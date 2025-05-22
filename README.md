
<!-- README.md is generated from README.Rmd. Please edit that file -->

## `apm`: Averaged Prediction Models

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/apm)](https://CRAN.R-project.org/package=apm)
<!-- badges: end -->

## Introduction

The `apm` package implements *Averaged Prediction Models (APM)*, a
Bayesian model averaging approach for controlled pre-post designs. These
designs compare differences over time between a group that becomes
exposed (treated group) and one that remains unexposed (comparison
group). With appropriate causal assumptions, they can identify the
causal effect of the exposure/treatment.

In APM, we specify a collection of models that predict untreated
outcomes. Our causal identifying assumption is that the model’s
prediction errors would be equal (in expectation) in the treated and
comparison groups in the absence of the exposure. This is a
generalization of familiar methods like Difference-in-Differences (DiD)
and Comparative Interrupted Time Series (CITS).

Because many models may be plausible for this prediction task, we
combine them using Bayesian model averaging. We weight each model by its
robustness to violations of the causal assumption.

## Installation

To install the package from CRAN, use

``` r
install.packages("apm")
```

To install the development version from GitHub, use:

``` r
# Install devtools if not already installed
install.packages("remotes")

remotes::install_github("tl2624/apm")
```

See `vignette("apm")` for details on using the package.
