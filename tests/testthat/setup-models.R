set.seed(123)
n <- 200

x1 <- rnorm(n)
x2 <- rnorm(n)

## logistic (binomial glm)
y <- rbinom(n, 1, plogis(0.8 * x1))
d <- data.frame(x1, x2, y)
model  <- glm(y ~ x1 + x2, data = d, family = binomial())
model0 <- glm(y ~ 1,       data = d, family = binomial())

## gaussian (lm) -- the classical-GLM benchmark
y_gauss <- 2 + 0.8 * x1 + rnorm(n)
d_gauss <- data.frame(x1, x2, y = y_gauss)
model_gauss  <- lm(y ~ x1 + x2, data = d_gauss)
model_gauss0 <- lm(y ~ 1,       data = d_gauss)

## poisson (glm)
y_pois <- rpois(n, exp(0.3 + 0.4 * x1))
d_pois <- data.frame(x1, x2, y = y_pois)
model_pois  <- glm(y ~ x1 + x2, data = d_pois, family = poisson())
model_pois0 <- glm(y ~ 1,       data = d_pois, family = poisson())

## multinomial (nnet::multinom); needs model = TRUE for r2()/eta2() to work
p_multinom <- plogis(0.8 * x1)
y_multinom <- apply(cbind(p_multinom, (1 - p_multinom) * 0.5, (1 - p_multinom) * 0.5), 1,
                     function(pr) sample(1:3, 1, prob = pr))
d_multinom <- data.frame(x1, x2, y = factor(y_multinom))
model_multinom  <- nnet::multinom(y ~ x1 + x2, data = d_multinom, trace = FALSE, model = TRUE)
model_multinom0 <- nnet::multinom(y ~ 1,       data = d_multinom, trace = FALSE, model = TRUE)

## ordinal (proportional-odds), fit two ways: MASS::polr and ordinal::clm
y_ord <- cut(0.8 * x1 + rnorm(n), breaks = c(-Inf, -0.5, 0.5, Inf), labels = 1:3)
d_ord <- data.frame(x1, x2, y = factor(y_ord, ordered = TRUE))

model_polr  <- MASS::polr(y ~ x1 + x2, data = d_ord, Hess = TRUE)
model_polr0 <- MASS::polr(y ~ 1,       data = d_ord, Hess = TRUE)

model_clm  <- ordinal::clm(y ~ x1 + x2, data = d_ord)
model_clm0 <- ordinal::clm(y ~ 1,       data = d_ord)
