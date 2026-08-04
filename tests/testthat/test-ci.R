test_that("ci=TRUE adds confidence intervals to the result structure", {
  expect_null(r2(model)$ci)
  expect_null(eta2(model)$ci)
  expect_null(eta2_partial(model)$ci)

  suppressWarnings({
    expect_equal(nrow(r2(model, ci=TRUE, ci_width=.90)$ci), 1)
    expect_equal(nrow(eta2(model, ci=TRUE, ci_width=.90)$ci), 2)
    expect_equal(nrow(eta2_partial(model, ci=TRUE, ci_width=.90)$ci), 2)
    expect_equal(nrow(eta2_partial(model_gauss, ci=TRUE, ci_width=.90)$ci), 2)
  })
})

test_that("partial confidence intervals use residual plus effect deviance", {
  res <- suppressWarnings(eta2_partial(model, ci=TRUE, ci_width=.90))
  expected <- ci_eta2p(
    eta2p=res$indices[1, "Eta2_p"],
    u=res$df[1],
    Dmx=res$D0 + res$D1[1],
    conf.level=.90
  )

  expect_equal(
    unname(unlist(res$ci[1, c("eta2p", "lower", "upper")])),
    unname(unlist(expected[, c("eta2p", "lower", "upper")]))
  )
  expect_error(r2(model, ci=TRUE, ci_width=1), "between 0 and 1")
})

test_that("quiet suppresses downstream warnings and messages", {
  expect_silent(r2(model, ci=TRUE, quiet=TRUE))
  expect_silent(eta2(model, ci=TRUE, quiet=TRUE))
  expect_silent(eta2_partial(model, ci=TRUE, quiet=TRUE))
})
