#library(gzlmpower)
set.seed(123)
n  <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
y  <- rbinom(n, 1, plogis(0.8 * x1))
d  <- data.frame(x1, x2, y)

model  <- glm(y ~   x2, data = d, family = binomial())
model0 <- glm(y ~ 1,       data = d, family = binomial())
# anova(model,model0)
eta2(model)
confint(a)


test_that("r2() matches the deviance-ratio formula for a logistic model", {
  res <- r2(model)

  k <- length(coef(model)) - length(coef(model0))
  expected_r2    <- 1 - deviance(model) / deviance(model0)
  expected_r2adj <- 1 - (deviance(model) + k) / deviance(model0)

  expect_equal(unname(res["r2"]),    expected_r2)
  expect_equal(unname(res["r2adj"]), expected_r2adj)
  expect_true(res["r2adj"] <= res["r2"])
  expect_true(res["r2"] >= 0 && res["r2"] <= 1)
})

test_that("eta2() matches the deviance-based formula for a logistic model", {
  res <- eta2(model)
  a   <- car::Anova(model, type = 3, test = "LR")

  expected_eta2 <- a[, "LR Chisq"] / deviance(model0)
  expected_eps2 <- pmax((a[, "LR Chisq"] - a$Df) / deviance(model0), 0)

  expect_equal(unname(res[, "Eta_squared"]),     unname(expected_eta2))
  expect_equal(unname(res[, "Epsilon_squared"]), unname(expected_eps2))
  expect_true(all(res[, "Epsilon_squared"] <= res[, "Eta_squared"]))
  expect_true(all(res >= 0 & res <= 1))
})
