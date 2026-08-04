test_that("r2() matches the deviance-ratio formula for an ordinal (polr) model", {
  res <- r2(model_polr)

  k <- length(coef(model_polr)) - length(coef(model_polr0))
  expected_r2    <- 1 - deviance(model_polr) / deviance(model_polr0)
  expected_r2adj <- 1 - (deviance(model_polr) + k) / deviance(model_polr0)

  expect_equal(res$indices$r2,    expected_r2)
  expect_equal(res$indices$r2adj, expected_r2adj)
  expect_true(res$indices$r2adj <= res$indices$r2)
  expect_true(res$indices$r2 >= 0 && res$indices$r2 <= 1)
})

test_that("eta2() matches the deviance-based formula for an ordinal (polr) model", {
  res <- eta2(model_polr)
  a   <- car::Anova(model_polr, type = 3)

  expected_eta2 <- a[, "LR Chisq"] / deviance(model_polr0)
  expected_eps2 <- pmax((a[, "LR Chisq"] - a$Df) / deviance(model_polr0), 0)

  expect_equal(unname(res$indices[, "Eta_squared"]),     unname(expected_eta2))
  expect_equal(unname(res$indices[, "Epsilon_squared"]), unname(expected_eps2))
  expect_true(all(res$indices[, "Epsilon_squared"] <= res$indices[, "Eta_squared"]))
  expect_true(all(res$indices >= 0 & res$indices <= 1))
})
