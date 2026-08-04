#library(gzlmpower)

n  <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
y  <- rbinom(n, 1, plogis(0.8 * x1))
d  <- data.frame(x1, x2, y)

lmodel  <- lm(y ~   x1+ x2, data = d)
res<-gzlmpower::r2(lmodel,test=T)
res
ss<-summary(lmodel)
ss$r.squared
gzlmpower::eta2(lmodel,ci=T,quiet=TRUE)
gzlmpower::r2(lmodel)

model  <- glm(y ~   x2, data = d, family = binomial())
model0 <- glm(y ~ 1,       data = d, family = binomial())
mod<-model
res<-gzlmpower::r2(mod,test=T)
res
res<-gzlmpower::r2(lmodel,test=T)
res
res$indices

res <- gzlmpower::r2(mod)
ci <- suppressWarnings(
  gzlmpower::ci_eta2(res$indices$r2,u=res$df,D0=res$D0)
)

test_that("r2() matches the deviance-ratio formula for a logistic model", {
  res <- r2(model)

  k <- length(coef(model)) - length(coef(model0))
  expected_r2    <- 1 - deviance(model) / deviance(model0)
  expected_r2adj <- 1 - (deviance(model) + k) / deviance(model0)

  expect_equal(unname(res$indices$r2),    expected_r2)
  expect_equal(unname(res$indices$r2adj), expected_r2adj)
  expect_true(res$indices$r2adj <= res$indices$r2)
  expect_true(res$indices$r2 >= 0 && res$indices$r2 <= 1)
})

test_that("eta2() matches the deviance-based formula for a logistic model", {
  res <- eta2(model)
  a   <- car::Anova(model, type = 3, test = "LR")

  expected_eta2 <- a[, "LR Chisq"] / deviance(model0)
  expected_eps2 <- pmax((a[, "LR Chisq"] - a$Df) / deviance(model0), 0)

  expect_equal(unname(res$indices[, "Eta_squared"]),
               unname(expected_eta2))
  expect_equal(unname(res$indices[, "Epsilon_squared"]),
               unname(expected_eps2))
  expect_true(all(res$indices[, "Epsilon_squared"] <=
                  res$indices[, "Eta_squared"]))
  expect_true(all(res$indices >= 0 & res$indices <= 1))
})
