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

test_that("eta2() matches the Chisq-based formula for an ordinal (clm) model", {
  res <- eta2(model_clm)
  a   <- car::Anova(model_clm, type = 3, test = "Chisq")

  expected_eta2 <- a[, "Chisq"] / dev_clm(model_clm0)
  expected_eps2 <- pmax((a[, "Chisq"] - a$Df) / dev_clm(model_clm0), 0)

  expect_equal(unname(res[, "Eta_squared"]),     unname(expected_eta2))
  expect_equal(unname(res[, "Epsilon_squared"]), unname(expected_eps2))
  expect_true(all(res[, "Epsilon_squared"] <= res[, "Eta_squared"]))
  expect_true(all(res >= 0 & res <= 1))
})
