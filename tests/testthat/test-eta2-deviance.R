test_that("eta2() stores null and effect deviances for GLMs", {
  res <- eta2(model)
  a <- car::Anova(model, type = 3, test = "LR")

  expect_equal(res$D0, deviance(model0))
  expect_equal(unname(res$D1), unname(a[, "LR Chisq"]))
  expect_equal(res$D, res$D1)
})

test_that("eta2() stores null and effect sums of squares for lm", {
  res <- eta2(model_gauss)
  a <- car::Anova(model_gauss, type = 3)
  w <- which(rownames(a) %in% c("(Intercept)", "Residuals"))
  effect_ss <- a$`Sum Sq`[-w]

  expect_equal(res$D0, deviance(model_gauss0))
  expect_equal(unname(res$D1), unname(effect_ss))
  expect_equal(res$D, res$D1)
})

test_that("eta2_partial() stores reduced and effect sums of squares for lm", {
  res <- eta2_partial(model_gauss)
  a <- car::Anova(model_gauss, type = 3)
  w <- which(rownames(a) %in% c("(Intercept)", "Residuals"))
  effect_ss <- a$`Sum Sq`[-w]

  expect_equal(res$D0, deviance(model_gauss))
  expect_equal(unname(res$D1), unname(effect_ss))
  expect_equal(res$D, res$D1)
  expect_true(is.null(res$anova))
  expect_true(!is.null(eta2_partial(model_gauss, test = TRUE)$anova))
})
