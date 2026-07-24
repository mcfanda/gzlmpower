# ordinal::clm has no deviance() method (it returns NULL), so the package falls back to
# -2*logLik(); the independent check here must do the same.
dev_clm <- function(object) -2 * as.numeric(logLik(object))

test_that("r2() matches the deviance-ratio formula for an ordinal (clm) model", {
  res <- r2(model_clm)

  k <- length(coef(model_clm)) - length(coef(model_clm0))
  expected_r2    <- 1 - dev_clm(model_clm) / dev_clm(model_clm0)
  expected_r2adj <- 1 - (dev_clm(model_clm) + k) / dev_clm(model_clm0)

  expect_equal(unname(res["r2"]),    expected_r2)
  expect_equal(unname(res["r2adj"]), expected_r2adj)
  expect_true(res["r2adj"] <= res["r2"])
  expect_true(res["r2"] >= 0 && res["r2"] <= 1)
})

test_that("eta2() matches a genuine per-term likelihood-ratio test for an ordinal (clm) model", {
  # car::Anova(clm_obj, type=3, test="Chisq") is NOT a likelihood-ratio test for `clm`
  # objects (Anova.clm just relabels the object and delegates to Anova.default, which is a
  # Wald test based on vcov()) -- so the independent reference here is built by manually
  # refitting the model with each term dropped in turn and comparing deviances directly,
  # rather than by calling car::Anova() or gzlmpower's own internal helper.
  res <- eta2(model_clm)

  model_clm_no_x1 <- ordinal::clm(y ~ x2, data = d_ord)
  model_clm_no_x2 <- ordinal::clm(y ~ x1, data = d_ord)

  chisq_x1 <- dev_clm(model_clm_no_x1) - dev_clm(model_clm)
  chisq_x2 <- dev_clm(model_clm_no_x2) - dev_clm(model_clm)
  df <- c(x1 = 1, x2 = 1)

  expected_eta2 <- c(chisq_x1, chisq_x2) / dev_clm(model_clm0)
  expected_eps2 <- pmax((c(chisq_x1, chisq_x2) - df) / dev_clm(model_clm0), 0)

  expect_equal(unname(res[, "Eta_squared"]),     unname(expected_eta2), tolerance = 1e-6)
  expect_equal(unname(res[, "Epsilon_squared"]), unname(expected_eps2), tolerance = 1e-6)
  expect_true(all(res[, "Epsilon_squared"] <= res[, "Eta_squared"]))
  expect_true(all(res >= 0 & res <= 1))
})

test_that("eta2_partial() matches a genuine per-term likelihood-ratio test for an ordinal (clm) model", {
  res <- eta2_partial(model_clm)

  model_clm_no_x1 <- ordinal::clm(y ~ x2, data = d_ord)
  model_clm_no_x2 <- ordinal::clm(y ~ x1, data = d_ord)

  chisq_x1 <- dev_clm(model_clm_no_x1) - dev_clm(model_clm)
  chisq_x2 <- dev_clm(model_clm_no_x2) - dev_clm(model_clm)
  chisq <- c(chisq_x1, chisq_x2)
  df <- c(x1 = 1, x2 = 1)
  k <- sum(df)

  devmx <- dev_clm(model_clm) + chisq   # deviance of the reduced (term-dropped) models
  expected_eta2 <- chisq / devmx
  expected_eps2 <- pmax((chisq - df) / (devmx + k - df), 0)

  expect_equal(unname(res[, "Eta_squared"]),     unname(expected_eta2), tolerance = 1e-6)
  expect_equal(unname(res[, "Epsilon_squared"]), unname(expected_eps2), tolerance = 1e-6)
  expect_true(all(res[, "Epsilon_squared"] <= res[, "Eta_squared"]))
  expect_true(all(res >= 0 & res <= 1))
})
