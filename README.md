
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FCBMA R package

Factor collapsing with Bayesian model averaging for regression models

Authors:

  - Sen Hu (<ethansen.hu@gmail.com>)
  - Adrian O’Hagan (<adrian.ohagan@ucd.ie>)
  - Thomas Brendan Murphy (<brendan.murphy@ucd.ie>)

## Description

The `FCBMA` package uses factor collapsing (FC) and Bayesian model
averaging (BMA) to find the optimal manners of combinations of
categorical levels within categorical variables (i.e. clulstering of
categorical levels) in linear or generalized linear regression models,
as introduced in Hu et al (2018)

## Installation

You can install the latest development version of `FCBMA` from
[GitHub](https://github.com/senhu/FCBMA):

``` r
install.packages("devtools")
devtools::install_github("senhu/FCBMA")
```

Then the package can be loaded with:

``` r
library(FCBMA)
```

## Example

This README file follows a package vignette format, and an example is
briefly demonstrated using the Swedish third party motor insurance
claims data in 1977 (available in this package), as illustrated in Hu et
al (2018). The data can be loaded via

``` r
data("sweden")
```

Details about the data set can be found in the data set documentation
within this package.

We start with the frequency aspect of claim modelling, by building a
baseline model where categorical variables are
unchanged:

``` r
freq <- glm(Claims ~ Zone + Bonus + Make + Kilometres, offset = log(Insured), data = sweden, family = "poisson")
```

Then the factors can be collpased individually, via

``` r
freq.kilo <- FCBMA(freq, varia.list = c("Kilometres"),
                   method = "complete", verbose = FALSE)
freq.make <- FCBMA(model = freq, varia.list = c("Make"),
                   method = "SA", verbose = FALSE)
```

or they can be collpased all together, via

``` r
freq.all <- FCBMA(freq, 
                  varia.list = c("Kilometres", "Zone", "Bonus", "Make"),
                  method = "GA", verbose = FALSE)
```

Depending on the size of the potential model space to search for the
opitmal combinations of collapsing, either greedy/complete search or
stochastic search with simulated annealing or genetic algorithm can be
used, by setting `method = "complete"`, `"SA"` or `"GA"`.

## Reference

Hu, S., O’Hagan, A., and Murphy, T. B. (2018) Motor insurance claim
modelling with factor collapsing and Bayesian model averaging. Stat, 7:
e180. doi: 10.1002/sta4.180.
