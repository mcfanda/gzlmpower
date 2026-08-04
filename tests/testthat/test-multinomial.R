test_that("r2() matches the deviance-ratio formula for a multinomial model", {
  res <- r2(model_multinom)

  k <- length(coef(model_multinom)) - length(coef(model_multinom0))
  expected_r2    <- 1 - deviance(model_multinom) / deviance(model_multinom0)
  expected_r2adj <- 1 - (deviance(model_multinom) + k) / deviance(model_multinom0)

  expect_equal(res$indices$r2,    expected_r2)
  expect_equal(res$indices$r2adj, expected_r2adj)
  expect_true(res$indices$r2adj <= res$indices$r2)
  expect_true(res$indices$r2 >= 0 && res$indices$r2 <= 1)
})

test_that("eta2() matches the deviance-based formula for a multinomial model", {
  res <- eta2(model_multinom)
  a   <- car::Anova(model_multinom, type = 3)

  expected_eta2 <- a[, "LR Chisq"] / deviance(model_multinom0)
  expected_eps2 <- pmax((a[, "LR Chisq"] - a$Df) / deviance(model_multinom0), 0)

  expect_equal(unname(res$indices[, "Eta_squared"]),     unname(expected_eta2))
  expect_equal(unname(res$indices[, "Epsilon_squared"]), unname(expected_eps2))
  expect_true(all(res$indices[, "Epsilon_squared"] <= res$indices[, "Eta_squared"]))
  expect_true(all(res$indices >= 0 & res$indices <= 1))
})

test_that("r2()/eta2() refuse a multinomial model fit without model = TRUE", {
  model_bad <- nnet::multinom(y ~ x1 + x2, data = d_multinom, trace = FALSE)

  expect_error(r2(model_bad),   "model=TRUE", fixed = TRUE)
  expect_error(eta2(model_bad), "model=TRUE", fixed = TRUE)
})
