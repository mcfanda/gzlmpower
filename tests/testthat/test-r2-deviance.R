test_that("r2() stores null and full-model deviances for GLMs", {
  res <- r2(model)

  expect_equal(res$D0, deviance(model0))
  expect_equal(res$D1, deviance(model))
  expect_equal(res$D, res$D1)
})

test_that("r2() stores null and full-model residual sums of squares for lm", {
  res <- r2(model_gauss)

  expect_equal(res$D0, deviance(model_gauss0))
  expect_equal(res$D1, deviance(model_gauss))
  expect_equal(res$D1, sum(residuals(model_gauss)^2))
  expect_equal(res$D, res$D1)
})

test_that("r2() reports a positive LR difference when test=TRUE", {
  res <- r2(model, test=TRUE)

  expect_equal(
    unname(res$anova$test),
    unname(deviance(model0) - deviance(model))
  )
  expect_gte(res$anova$test, 0)
})
