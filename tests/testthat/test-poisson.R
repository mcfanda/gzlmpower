test_that("r2() matches the deviance-ratio formula for a poisson model", {
  res <- r2(model_pois)

  k <- length(coef(model_pois)) - length(coef(model_pois0))
  expected_r2    <- 1 - deviance(model_pois) / deviance(model_pois0)
  expected_r2adj <- 1 - (deviance(model_pois) + k) / deviance(model_pois0)

  expect_equal(unname(res["r2"]),    expected_r2)
  expect_equal(unname(res["r2adj"]), expected_r2adj)
  expect_true(res["r2adj"] <= res["r2"])
  expect_true(res["r2"] >= 0 && res["r2"] <= 1)
})

test_that("eta2() matches the deviance-based formula for a poisson model", {
  res <- eta2(model_pois)
  a   <- car::Anova(model_pois, type = 3, test = "LR")

  expected_eta2 <- a[, "LR Chisq"] / deviance(model_pois0)
  expected_eps2 <- pmax((a[, "LR Chisq"] - a$Df) / deviance(model_pois0), 0)

  expect_equal(unname(res[, "Eta_squared"]),     unname(expected_eta2))
  expect_equal(unname(res[, "Epsilon_squared"]), unname(expected_eps2))
  expect_true(all(res[, "Epsilon_squared"] <= res[, "Eta_squared"]))
  expect_true(all(res >= 0 & res <= 1))
})
