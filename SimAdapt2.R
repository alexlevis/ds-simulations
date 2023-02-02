library(tidyverse)
library(mvtnorm)
library(xtable)

dat <- read.csv("durable_missing_cbl_only.csv",
                stringsAsFactors = T)
dat$A <- as.numeric(dat$final_decision_merge) - 1
dat$L <- 1*(dat$GENDER == "M")

## try and understand outcome distribution so simulation makes sense
summary(dat$ptwc)
mean(dat$ptwc[dat$A == 0]) # -0.236
mean(dat$ptwc[dat$A == 1]) # -0.171
hist(dat$ptwc)

n <- nrow(dat)
### setting up simulation parameters
## A | L
p1 <- mean(dat$A[dat$L == 1]) + 0.06
p0 <- mean(dat$A[dat$L == 0]) - 0.06

## R | L, A
delta0 <- qlogis(0.2) # marginal probability of missing outcome approx 20%
deltaL <- 0.09
deltaA <- -0.05
deltaLA <- -0.35

## Y | L, A, R
outcome.true <- lm(ptwc ~ L + A, data = dat)
summary(outcome.true)
beta0 <- coef(outcome.true)[1]
betaA <- coef(outcome.true)[3]
betaL <- coef(outcome.true)[2]
betaRA <- 0.000 ### PLAY WITH THIS TO MODIFY MNAR VIOLATION
sig.star <- sigma(outcome.true)

## S | L, A, R = 0, try to make it so that get a good number of each group
zeta0 <- qlogis(0.10) # marginal probability of double sampling given R == 0 approx 15%
zetaL <- 0.4
zetaA <- 0.3
zetaLA <- 0.25

### DETERMINE TRUTH
N <- 20000
## dist of L sampled from the empirical distribution
dat.true <- cbind.data.frame(L = sample(dat$L, N, replace = T))

## dist of A | L
dat.true$A <- rbinom(N, 1, ifelse(dat.true$L == 1, p1, p0))

## dist of R | L, A
dat.true$R <- rbinom(N, 1, plogis(delta0 + deltaA * dat.true$A +
                                    deltaL * dat.true$L +
                                    deltaLA * dat.true$L * dat.true$A))

chi.0 <- mean(beta0 + betaL* dat.true$L)
chi.1 <- mean((beta0 + betaA + betaL * dat.true$L + betaRA) * 
                plogis(delta0 + deltaA  + deltaL * dat.true$L + deltaLA * dat.true$L)) +
  mean((beta0 + betaA + betaL * dat.true$L) * 
         (1 - plogis(delta0 + deltaA  + deltaL * dat.true$L + deltaLA * dat.true$L)))

# simulation time!
set.seed(145)
R <- 5000 ## number of simulations
res <- matrix(NA, nrow = R, ncol = 4)
res.sig <- matrix(NA, nrow = R, ncol = 5)
choice.adapt <- choice.naive <- rep(NA, R)
for (r in 1:R) {
  ## dist of L sampled from the empirical distribution
  dat.sim <- cbind.data.frame(L = sample(dat$L, n, replace = T))
  
  ## dist of A | L
  dat.sim$A <- rbinom(n, 1, ifelse(dat.sim$L == 1, p1, p0))
  
  ## dist of R | L, A
  dat.sim$R <- rbinom(n, 1, plogis(delta0 + deltaA * dat.sim$A +
                                     deltaL * dat.sim$L +
                                     deltaLA * dat.sim$L * dat.sim$A))
  
  ## dist of Y | L, A, R
  dat.sim$Y <- 
    rnorm(n, 
          beta0 + betaA * dat.sim$A + betaL* dat.sim$L +
            betaRA * dat.sim$R * dat.sim$A, sig.star)
  
  ## dist of S | L, A, R = 0
  eta.1 <- plogis(zeta0 + zetaA * dat.sim$A +
                    zetaL * dat.sim$L +
                    zetaLA * dat.sim$L * dat.sim$A)
  dat.sim$S <- rbinom(n, 1, ifelse(dat.sim$R == 1, 0, eta.1))
  
  ## fitting models
  pi.hat <- glm(A ~ L, data = dat.sim)
  gamma.hat <- glm(R ~ A * L, data = dat.sim)
  muR.hat <- lm(Y ~ A + L, data = dat.sim[dat.sim$R == 1,])
  muS.hat <- lm(Y ~ A + L, data = dat.sim[dat.sim$S == 1,])
  
  ## for MAR.eff
  muRS.hat <- lm(Y ~ A + L, data = dat.sim[dat.sim$R == 1 | dat.sim$S == 1,])
  
  dat.a1 <- cbind.data.frame(L = dat.sim$L, A = 1)
  dat.a0 <- cbind.data.frame(L = dat.sim$L, A = 0)
  ## estimators
  mu.1.a1 <- predict(muR.hat, newdata = dat.a1)
  mu.0.a1 <- predict(muS.hat, newdata = dat.a1)
  mu.1.a0 <- predict(muR.hat, newdata = dat.a0)
  mu.0.a0 <- predict(muS.hat, newdata = dat.a0)
  pi.1 <- predict(pi.hat, newdata = dat.sim)
  gamma.a1 <- predict(gamma.hat, newdata = dat.a1)
  gamma.a0 <- predict(gamma.hat, newdata = dat.a0)
  mu.a1 <- mu.1.a1 * gamma.a1 + mu.0.a1 * (1 - gamma.a1)
  mu.a0 <- mu.1.a0 * gamma.a0 + mu.0.a0 * (1 - gamma.a0)
  
  # for MAR.eff
  mu.e.a1 <- predict(muRS.hat, newdata = dat.a1)
  mu.e.a0 <- predict(muRS.hat, newdata = dat.a0)
  
  IF.0 <- mu.a0 + 
    ifelse(dat.sim$A == 0, 
           ((dat.sim$R + dat.sim$S / eta.1) * dat.sim$Y - mu.a0 +
              (1 - dat.sim$R) * (1 - dat.sim$S / eta.1) * mu.0.a0) / (1 - pi.1),
           0)
  IF.1 <- mu.a1 + 
    ifelse(dat.sim$A == 1, 
           ((dat.sim$R + dat.sim$S / eta.1) * dat.sim$Y - mu.a1 +
              (1 - dat.sim$R) * (1 - dat.sim$S / eta.1) * mu.0.a1) / pi.1,
           0)
  IF.MAR.eff.0 <- mu.e.a0 + (dat.sim$R + dat.sim$S) * (1 - dat.sim$A) *
    (dat.sim$Y - mu.e.a0) / ((1 - pi.1) * (gamma.a0 + (1 - gamma.a0) * eta.1))
  IF.MAR.eff.1 <- mu.e.a1 + (dat.sim$R + dat.sim$S) * dat.sim$A *
    (dat.sim$Y - mu.e.a1) / (pi.1 * (gamma.a1 + (1 - gamma.a1) * eta.1))
  
  ## for computing confidence intervals and Dominik Rothenhausler's criterion
  sigma.0.2 <- mean((IF.1 - IF.0 - mean(IF.1 - IF.0))^2)
  sigma.1.2 <- mean((IF.MAR.eff.1 - IF.MAR.eff.0 - 
                       mean(IF.MAR.eff.1 - IF.MAR.eff.0))^2)
  tau.1.2 <- mean((IF.1 - IF.0 - IF.MAR.eff.1 + IF.MAR.eff.0 - 
                     mean(IF.1 - IF.0 - IF.MAR.eff.1 + IF.MAR.eff.0))^2)
  sigma.cov01 <- mean((IF.1 - IF.0 - mean(IF.1 - IF.0)) *
                        (IF.MAR.eff.1 - IF.MAR.eff.0 - 
                           mean(IF.MAR.eff.1 - IF.MAR.eff.0)))
  ## compute Dominik Rothenhausler criterion
  # R.hat <- c(sigma.0.2 / n, (sigma.1.2 - tau.1.2)/n + 
  #              (mean((IF.MAR.eff.1 - IF.MAR.eff.0) - (IF.1 - IF.0)))^2)
  R.mod <- c(sigma.0.2 / n, sigma.1.2/n + 
               max(mean(((IF.MAR.eff.1 - IF.MAR.eff.0) - (IF.1 - IF.0)))^2 - 
                     tau.1.2/n, 0))
  
  ## compute Wald test statistic and ensuing naive estimator
  Wald.test <- 
    sqrt(n) * 
    (mean(IF.1) - mean(IF.0) - (mean(IF.MAR.eff.1) - mean(IF.MAR.eff.0))) /
    sqrt(sigma.0.2 + sigma.1.2 - 2 * sigma.cov01)
  
  choice.adapt[r] <- 1*(R.mod[1] < R.mod[2])
  choice.naive[r] <- 1*(abs(Wald.test) > qnorm(0.975))
  
  res[r,] <- c(mean(IF.1) - mean(IF.0),
               mean(IF.MAR.eff.1) - mean(IF.MAR.eff.0),
               ifelse(R.mod[1] < R.mod[2], mean(IF.1) - mean(IF.0),
                      mean(IF.MAR.eff.1) - mean(IF.MAR.eff.0)),
               ifelse(abs(Wald.test) > qnorm(0.975), mean(IF.1) - mean(IF.0),
                      mean(IF.MAR.eff.1) - mean(IF.MAR.eff.0)))
  res.sig[r,] <- c(sigma.0.2, sigma.1.2, tau.1.2, sigma.cov01, 
                   ifelse(abs(Wald.test) > qnorm(0.975), sigma.0.2, sigma.1.2))
}

covrg.adapt.hits <- rep(NA, R)
length.adapt.hits <- rep(NA, R)

## compute CIs for adaptive estimator
set.seed(264)
M <- 5000
beta.grid <- seq(-2, 2, by = 0.02)
for (r in 1:R) {
  Sigma <- rbind(c(res.sig[r, 1], res.sig[r, 4]),
                 c(res.sig[r, 4], res.sig[r, 2]))
  Z <- rmvnorm(M, sqrt(n) * res[r, c(1, 2)], Sigma)
  ## compute indicator of the event A_0^{est}
  A.g <- 1 * (res.sig[r, 1] < pmax((Z[,2] - Z[,1])^2, 0) + res.sig[r, 2])
  b.beta <- sapply(beta.grid, function(x) {
    mean((Z[,1] - sqrt(n) * res[r, 1] <= x) & A.g == 1) +
      mean((Z[,2] - sqrt(n) * res[r, 1] <= x) & A.g == 0)
  }, simplify = 0)
  b.025 <- beta.grid[which.min(abs(b.beta - 0.025))]
  b.975 <- beta.grid[which.min(abs(b.beta - 0.975))]
  covrg.adapt.hits[r] <- 1 * ((chi.1 - chi.0 >= res[r, 3] - b.975 / sqrt(n)) &
                                (chi.1 - chi.0 <= res[r, 3] - b.025 / sqrt(n)))
  length.adapt.hits[r] <- (b.975 - b.025) / sqrt(n)
}

covrg.NP <- mean(abs(res[,1] - (chi.1 - chi.0)) <= 
                   qnorm(0.975) * sqrt(res.sig[,1] / n))
length.NP <- mean(2 * qnorm(0.975) * sqrt(res.sig[,1] / n))
covrg.MAR.eff <- mean(abs(res[,2] - (chi.1 - chi.0)) <= 
                        qnorm(0.975) * sqrt(res.sig[,2] / n))
length.MAR.eff <- mean(2 * qnorm(0.975) * sqrt(res.sig[,2] / n))
covrg.naive <- mean(abs(res[,4] - (chi.1 - chi.0)) <= 
                      qnorm(0.975) * sqrt(res.sig[,5] / n))
length.naive <- mean(2 * qnorm(0.975) * sqrt(res.sig[,5] / n))
covrg.adapt <- mean(covrg.adapt.hits)
length.adapt <- mean(length.adapt.hits)

# save(res, res.sig, choice.adapt, choice.naive,
#      covrg.NP, covrg.MAR.eff, covrg.naive, covrg.adapt,
#      length.NP, length.MAR.eff, length.naive, length.adapt,
#      file = "RES032.Rdata")
load("RES032.RData")
RESULTS032 <- c(mean(res[,1] - (chi.1 - chi.0)), var(res[,1]),
                mean(res[,2] - (chi.1 - chi.0)), var(res[,2]),
                mean(res[,3] - (chi.1 - chi.0)), var(res[,3]),
                mean(res[,4] - (chi.1 - chi.0)), var(res[,4]),
                covrg.NP, covrg.MAR.eff, covrg.adapt, covrg.naive,
                length.NP, length.MAR.eff, length.adapt, length.naive,
                mean(choice.adapt), mean(choice.naive))
round(RESULTS032[17:18], 2) # proportion of times choosing base estimator
print032 <- cbind(round(100 * RESULTS032[c(1,3,5,7)] / (chi.1 - chi.0), 2), # % bias
                  round(RESULTS032[c(2,4,6,8)] / RESULTS032[2], 2), # relative variance
                  round((RESULTS032[c(1,3,5,7)]^2 + RESULTS032[c(2,4,6,8)]) / 
                          (RESULTS032[1]^2 + RESULTS032[2]), 2), # relative MSE
                  round(100 * RESULTS032[9:12], 1), # coverage
                  round(RESULTS032[13:16], 3)) # length
xtable(print032, digits = 3)

# save(res, res.sig, choice.adapt, choice.naive,
#      covrg.NP, covrg.MAR.eff, covrg.naive, covrg.adapt,
#      length.NP, length.MAR.eff, length.naive, length.adapt,
#      file = "RES016.Rdata")
load("RES016.RData")
RESULTS016 <- c(mean(res[,1] - (chi.1 - chi.0)), var(res[,1]),
                mean(res[,2] - (chi.1 - chi.0)), var(res[,2]),
                mean(res[,3] - (chi.1 - chi.0)), var(res[,3]),
                mean(res[,4] - (chi.1 - chi.0)), var(res[,4]),
                covrg.NP, covrg.MAR.eff, covrg.adapt, covrg.naive,
                length.NP, length.MAR.eff, length.adapt, length.naive,
                mean(choice.adapt), mean(choice.naive))
round(RESULTS016[17:18], 2) # proportion of times choosing base estimator
print016 <- cbind(round(100 * RESULTS016[c(1,3,5,7)] / (chi.1 - chi.0), 2), # % bias
                  round(RESULTS016[c(2,4,6,8)] / RESULTS016[2], 2), # relative variance
                  round((RESULTS016[c(1,3,5,7)]^2 + RESULTS016[c(2,4,6,8)]) / 
                          (RESULTS016[1]^2 + RESULTS016[2]), 2), # relative MSE
                  round(100 * RESULTS016[9:12], 1), # coverage
                  round(RESULTS016[13:16], 3)) # length
xtable(print016, digits = 3)

# save(res, res.sig, choice.adapt, choice.naive,
#      covrg.NP, covrg.MAR.eff, covrg.naive, covrg.adapt,
#      length.NP, length.MAR.eff, length.naive, length.adapt,
#      file = "RES000.Rdata")
load("RES000.RData")
RESULTS000 <- c(mean(res[,1] - (chi.1 - chi.0)), var(res[,1]),
                mean(res[,2] - (chi.1 - chi.0)), var(res[,2]),
                mean(res[,3] - (chi.1 - chi.0)), var(res[,3]),
                mean(res[,4] - (chi.1 - chi.0)), var(res[,4]),
                covrg.NP, covrg.MAR.eff, covrg.adapt, covrg.naive,
                length.NP, length.MAR.eff, length.adapt, length.naive,
                mean(choice.adapt), mean(choice.naive))
round(RESULTS000[17:18], 2) # proportion of times choosing base estimator
print000 <- cbind(round(100 * RESULTS000[c(1,3,5,7)] / (chi.1 - chi.0), 2), # % bias
                  round(RESULTS000[c(2,4,6,8)] / RESULTS000[2], 2), # relative variance
                  round((RESULTS000[c(1,3,5,7)]^2 + RESULTS000[c(2,4,6,8)]) / 
                          (RESULTS000[1]^2 + RESULTS000[2]), 2), # relative MSE
                  round(100 * RESULTS000[9:12], 1), # coverage
                  round(RESULTS000[13:16], 3)) # length
xtable(print000, digits = 3)

# RES <- rbind(RESULTS000, RESULTS016, RESULTS032)
# colnames(RES) <- c("MAR", "Low MNAR", "High MNAR")
# rownames(RES) <- c("NP", "MAR eff", "Adaptive", "Naive Adaptive")
# RES
