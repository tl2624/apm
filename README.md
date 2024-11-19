
<!-- README.md is generated from README.Rmd. Please edit that file -->

# eepd

<!-- badges: start -->
<!-- badges: end -->

``` r
library(eepd)
data("ptpdata")
```

## Supplying models

We can specify the models to test using `eepd_mod()`. This create a full
cross of all supplied arguments, which include model formula, families,
whether the outcome is logged or not, whether fixed effects are included
or not, whether the outcome should be a difference, and whether outcome
lags should appear as predictors. Below, we create a cross of 9 models.

``` r
models <- eepd_mod(deaths ~ 1,
                   family = list("gaussian", "quasipoisson"),
                   log = c(TRUE, FALSE),
                   lag = 0, diff_k = 0, 
                   time_trend = 0:2)

models
#> - Model 1: baseline mean (log, Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: none
#> log outcome: yes
#> time trend: none
#> unit fixed effects: no
#> 
#> - Model 2: baseline mean (Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: none
#> log outcome: no
#> time trend: none
#> unit fixed effects: no
#> 
#> - Model 3: baseline mean (Quasipoisson)
#> deaths ~ 1
#> family: quasipoisson(link = 'log')
#> outcome lag: none
#> outcome diff: none
#> log outcome: no
#> time trend: none
#> unit fixed effects: no
#> 
#> - Model 4: linear trend (log, Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: none
#> log outcome: yes
#> time trend: linear
#> unit fixed effects: no
#> 
#> - Model 5: linear trend (Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: none
#> log outcome: no
#> time trend: linear
#> unit fixed effects: no
#> 
#> - Model 6: linear trend (Quasipoisson)
#> deaths ~ 1
#> family: quasipoisson(link = 'log')
#> outcome lag: none
#> outcome diff: none
#> log outcome: no
#> time trend: linear
#> unit fixed effects: no
#> 
#> - Model 7: quadratic trend (log, Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: none
#> log outcome: yes
#> time trend: quadratic
#> unit fixed effects: no
#> 
#> - Model 8: quadratic trend (Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: none
#> log outcome: no
#> time trend: quadratic
#> unit fixed effects: no
#> 
#> - Model 9: quadratic trend (Quasipoisson)
#> deaths ~ 1
#> family: quasipoisson(link = 'log')
#> outcome lag: none
#> outcome diff: none
#> log outcome: no
#> time trend: quadratic
#> unit fixed effects: no
```

Normally, this cross would yield 12 = 3 (formulas) x 2 (families) x 2
(log T/F), but by default any models with non-linear links and
`log = TRUE` are removed, leaving 9 models. If we want to manually add
other models, we can so by creating a new models object and appending it
to the current one.

``` r
models2 <- eepd_mod(list(deaths ~ 1),
                    diff_k = 1)

models <- c(models, models2)

models
#> - Model 1: baseline mean (log, Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: none
#> log outcome: yes
#> time trend: none
#> unit fixed effects: no
#> 
#> - Model 2: baseline mean (Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: none
#> log outcome: no
#> time trend: none
#> unit fixed effects: no
#> 
#> - Model 3: baseline mean (Quasipoisson)
#> deaths ~ 1
#> family: quasipoisson(link = 'log')
#> outcome lag: none
#> outcome diff: none
#> log outcome: no
#> time trend: none
#> unit fixed effects: no
#> 
#> - Model 4: linear trend (log, Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: none
#> log outcome: yes
#> time trend: linear
#> unit fixed effects: no
#> 
#> - Model 5: linear trend (Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: none
#> log outcome: no
#> time trend: linear
#> unit fixed effects: no
#> 
#> - Model 6: linear trend (Quasipoisson)
#> deaths ~ 1
#> family: quasipoisson(link = 'log')
#> outcome lag: none
#> outcome diff: none
#> log outcome: no
#> time trend: linear
#> unit fixed effects: no
#> 
#> - Model 7: quadratic trend (log, Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: none
#> log outcome: yes
#> time trend: quadratic
#> unit fixed effects: no
#> 
#> - Model 8: quadratic trend (Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: none
#> log outcome: no
#> time trend: quadratic
#> unit fixed effects: no
#> 
#> - Model 9: quadratic trend (Quasipoisson)
#> deaths ~ 1
#> family: quasipoisson(link = 'log')
#> outcome lag: none
#> outcome diff: none
#> log outcome: no
#> time trend: quadratic
#> unit fixed effects: no
#> 
#> - Model 10: baseline mean (1st diff, Gaussian)
#> deaths ~ 1
#> family: gaussian(link = 'identity')
#> outcome lag: none
#> outcome diff: 1
#> log outcome: no
#> time trend: none
#> unit fixed effects: no
```

This leaves us with 10 models.

## Fitting the models

Next we fit all 10 models to the data. We do so once for each validation
time to compute the average prediction error that will be used to select
the optimal model. All models are fit simultaneously so the simulation
can use the full joint distribution of model parameter estimates. For
each validation time, each model is fit using a dataset that contains
data points prior to that time.

We use `eepd_fit()` to fit the models, and calculate the prediction
errors and BMA weights.

``` r
fits <- eepd_pre(models,
                 data = ptpdata,
                 group_var = "group",
                 time_var = "year",
                 unit_var = "state",
                 val_times = 2004:2007)
#> Fitting models... Done.
#> Simulating to compute BMA weights...
#> Done.

fits
#> An `eepd_pre_fits` object
#> 
#>  - grouping variable: group
#>  - unit variable: state
#>  - time variable: year
#>    - validation times: 2004, 2005, 2006, 2007
#>  - number of models compared: 10
#>  - number of simulation iterations: 1000
#> 
#> Use `summary()` or `plot()` to examine prediction errors and BMA weights.
```

## Computing the ATT

We compute the ATT using `eepd_est()`, which uses bootstrapping to
compute model uncertainty due to sampling along with uncertainty due to
model selection.

``` r
est <- eepd_est(fits,
                post_time = 2008,
                M = 1,
                R = 50)

est
#> An `eepd_est` object
#> 
#>  - grouping variable: group
#>  - unit variable: state
#>  - time variable: year
#>    - validation times: 
#>    - post-treatment time: 2008
#>  - sensitivity parameter (M): 1
#>  - bootstrap replications: 50
#> 
#> Use `summary()` or `plot()` to examine estimates and uncertainty bounds.

summary(est)
#>       Estimate Std. Error  CI low CI high z_value Pr(>|z|)
#> ATT      62.00      40.45  -17.28  141.28   1.533    0.125
#> M = 1        .          . -116.69  181.35       .        .
```
