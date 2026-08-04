# gzlmpower

Model-based effect sizes and likelihood-ratio power analysis for generalized linear models.

## Overview

`gzlmpower` is a small R package for two closely related jobs:

- estimating effect sizes from fitted models;
- using those effect sizes in likelihood-ratio test power analysis.

The package focuses on generalized linear modeling workflows where researchers often want a compact summary such as:

- model-level `R^2`;
- term-level `eta^2` and `epsilon^2`;
- partial `eta^2` and partial `epsilon^2`;
- required sample size, power, alpha, or effect size for a likelihood-ratio test.

For non-Gaussian models, these quantities are based on deviance and likelihood-ratio tests rather than ordinary least-squares sums of squares.

## Installation

`gzlmpower` is not on CRAN. Install it from GitHub:

```r
pak::pak("mcfanda/gzlmpower")
# or
remotes::install_github("mcfanda/gzlmpower")
```

Some model classes used in the examples below come from suggested packages:

- `MASS` for `polr()`
- `nnet` for `multinom()`
- `ordinal` for `clm()`
- `MBESS` for confidence intervals

## What It Provides

### `r2()`

Returns model-level `R^2` and adjusted `R^2`.

For GLMs and related likelihood-based models, `gzlmpower` computes `R^2` from the intercept-only model and the fitted model:

```text
R^2 = 1 - D_model / D_null
```

### `eta2()`

Returns term-level `Eta_squared` and `Epsilon_squared`.

- For `lm`, these are based on Type III sums of squares.
- For `glm`, `MASS::polr`, and `nnet::multinom`, they are based on Type III likelihood-ratio tests.
- For `ordinal::clm`, the package refits reduced models term by term to obtain genuine likelihood-ratio tests, because the usual single-model ANOVA output is Wald-based.

### `eta2_partial()`

Returns partial `Eta_squared` and partial `Epsilon_squared`, using the reduced model that omits each term in turn.

### `power.lrt()`

Solves one missing quantity among:

- effect size (`es`)
- sample size (`N`)
- degrees of freedom (`df`)
- significance level (`sig.level`)
- power

The function is intended for likelihood-ratio or chi-square style tests for categorical outcomes. You must leave exactly one of those arguments as `NULL`.

## Supported Model Classes

The package is tested on:

- `lm`
- `glm` including logistic and Poisson models
- `MASS::polr`
- `ordinal::clm`
- `nnet::multinom`

For `multinom()` models, fit with `model = TRUE`, otherwise `r2()`, `eta2()`, and `eta2_partial()` cannot reconstruct the null model.

## Quick Start

The package ships with a small simulated dataset, `manymodels`, designed for trying different model families.

```r
library(gzlmpower)
data(manymodels)
```

### Logistic regression

```r
fit_binom <- glm(ybin ~ x + z, family = binomial(), data = manymodels)

r2(fit_binom)
eta2(fit_binom)
eta2_partial(fit_binom)
```

This gives a model-level `R^2` plus effect-size estimates for each predictor based on likelihood-ratio tests.

### Multinomial regression

```r
manymodels$ycat <- factor(manymodels$ycat)

fit_multi <- nnet::multinom(
  ycat ~ x + z,
  data = manymodels,
  trace = FALSE,
  model = TRUE
)

r2(fit_multi)
eta2(fit_multi)
```

### Ordinal regression

```r
manymodels$yord <- ordered(manymodels$yord)

fit_ord <- ordinal::clm(yord ~ x + z, data = manymodels)

r2(fit_ord)
eta2(fit_ord)
eta2_partial(fit_ord)
```

If you prefer `MASS::polr()`, that works too:

```r
fit_polr <- MASS::polr(yord ~ x + z, data = manymodels, Hess = TRUE)

r2(fit_polr)
eta2(fit_polr)
```

## Power Analysis

`power.lrt()` uses an expected effect size and the expected outcome distribution to solve for the missing quantity.

For a binary outcome, `prob` can be a single proportion:

```r
power.lrt(
  es = 0.10,
  prob = 0.50,
  df = 1,
  power = 0.80
)
```

The call above returns the required sample size `N`.

For outcomes with more than two categories, pass the full vector of expected proportions:

```r
power.lrt(
  es = 0.08,
  prob = c(0.20, 0.50, 0.30),
  df = 2,
  N = 180,
  sig.level = 0.05
)
```

## Confidence Intervals

Set `ci = TRUE` to add noncentral-chi-square confidence intervals to the result. Use `ci_width` to control the confidence level:

```r
fit_results <- eta2(fit_binom, ci = TRUE, ci_width = 0.95)
fit_results$ci
print(fit_results)
```

The intervals are returned in `$ci` and are printed beside the primary effect-size index. For `eta2_partial()`, the interval calculation uses the residual deviance or residual sum of squares plus the corresponding effect contribution.

## Notes

- `eta2()` and `eta2_partial()` rely on Type III tests through `car::Anova()` when that is appropriate for the model class.
- `ordinal::clm` is handled separately to obtain likelihood-ratio tests instead of Wald tests.
- `nnet::multinom()` models must be fit with `model = TRUE`.
- For binary outcomes, `prob = 0.5` is interpreted as `c(0.5, 0.5)`.

## License

GPL-3.
