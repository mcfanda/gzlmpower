test_that("r2() matches base R's r.squared for a gaussian (lm) model", {
  res <- r2(model_gauss)
  ss  <- summary(model_gauss)

  expect_equal(res$indices$r2,    ss$r.squared)
  expect_equal(res$indices$r2adj, ss$adj.r.squared)
  expect_true(res$indices$r2adj <= res$indices$r2)
  expect_true(res$indices$r2 >= 0 && res$indices$r2 <= 1)
})

test_that("eta2() matches the Type III sum-of-squares formula for a gaussian (lm) model", {
  res <- eta2(model_gauss)
  a   <- car::Anova(model_gauss, type = 3)
  w   <- which(rownames(a) %in% c("(Intercept)", "Residuals"))

  sse0 <- deviance(model_gauss0)
  msem <- deviance(model_gauss) / model_gauss$df.residual

  expected_eta2 <- a$`Sum Sq`[-w] / sse0
  expected_eps2 <- pmax((a$`Sum Sq`[-w] - a$Df[-w] * msem) / sse0, 0)

  expect_equal(unname(res$indices[, "Eta_squared"]),     unname(expected_eta2))
  expect_equal(unname(res$indices[, "Epsilon_squared"]), unname(expected_eps2))
  expect_true(all(res$indices[, "Epsilon_squared"] <= res$indices[, "Eta_squared"]))
  expect_true(all(res$indices >= 0 & res$indices <= 1))
})
