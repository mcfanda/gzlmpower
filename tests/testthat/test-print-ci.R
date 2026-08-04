test_that("print.nwa places confidence limits beside the primary index", {
  object <- structure(
    list(
      indices=data.frame(r2=0.20, r2adj=0.18, row.names="Model"),
      ci=data.frame(r2=0.20, lower=0.10, upper=0.30, row.names="Model"),
      anova=NULL
    ),
    class="nwa"
  )

  output <- capture.output(print(object))
  expect_true(any(grepl("CI_lower", output, fixed=TRUE)))
  expect_true(any(grepl("CI_upper", output, fixed=TRUE)))
})
