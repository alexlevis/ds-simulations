library(tidyverse)

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
betaRA <- 0.016 ### PLAY WITH THIS TO MODIFY MNAR VIOLATION
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

set.seed(145)
R <- 5000 ## number of simulations
res <- matrix(NA, nrow = R, ncol = 16)
colnames(res) <- c("IF.0", "IF.1", 
                   "IF.0pw", "IF.1pw",
                   "IF.0mw", "IF.1mw",
                   "IF.0bw", "IF.1bw",
                   "OR.MAR.0", "OR.MAR.1",
                   "IPW.MAR.0", "IPW.MAR.1",
                   "IF.MAR.0", "IF.MAR.1",
                   "IF.MAR.eff.0", "IF.MAR.eff.1")
tau.1.2 <- sigma.1.2 <- sigma.0.2 <- rep(NA, R)
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
  pi.wrong <- glm(1 - A ~ L, data = dat.sim)
  gamma.hat <- glm(R ~ A * L, data = dat.sim)
  gamma.wrong <- glm(R ~ L, data = dat.sim)
  muR.hat <- lm(Y ~ A + L, data = dat.sim[dat.sim$R == 1,])
  muS.hat <- lm(Y ~ A + L, data = dat.sim[dat.sim$S == 1,])
  muR.wrong <- lm(Y ~ A, data = dat.sim[dat.sim$R == 1,])
  muS.wrong <- lm(Y ~ A, data = dat.sim[dat.sim$S == 1,])
  
  ## for MAR.eff
  muRS.hat <- lm(Y ~ A + L, data = dat.sim[dat.sim$R == 1 | dat.sim$S == 1,])
  
  dat.a1 <- cbind.data.frame(L = dat.sim$L, A = 1)
  dat.a0 <- cbind.data.frame(L = dat.sim$L, A = 0)
  ## estimators
  mu.1.a1 <- predict(muR.hat, newdata = dat.a1)
  mu.0.a1 <- predict(muS.hat, newdata = dat.a1)
  mu.1.a0 <- predict(muR.hat, newdata = dat.a0)
  mu.0.a0 <- predict(muS.hat, newdata = dat.a0)
  mu.1.a1w <- predict(muR.wrong, newdata = dat.a1)
  mu.0.a1w <- predict(muS.wrong, newdata = dat.a1)
  mu.1.a0w <- predict(muR.wrong, newdata = dat.a0)
  mu.0.a0w <- predict(muS.wrong, newdata = dat.a0)
  pi.1 <- predict(pi.hat, newdata = dat.sim)
  pi.1w <- predict(pi.wrong, newdata = dat.sim)
  gamma.a1 <- predict(gamma.hat, newdata = dat.a1)
  gamma.a0 <- predict(gamma.hat, newdata = dat.a0)
  gamma.a1w <- predict(gamma.wrong, newdata = dat.a1)
  gamma.a0w <- predict(gamma.wrong, newdata = dat.a0)
  mu.a1 <- mu.1.a1 * gamma.a1 + mu.0.a1 * (1 - gamma.a1)
  mu.a0 <- mu.1.a0 * gamma.a0 + mu.0.a0 * (1 - gamma.a0)
  mu.a1w <- mu.1.a1w * gamma.a1 + mu.0.a1w * (1 - gamma.a1)
  mu.a0w <- mu.1.a0w * gamma.a0 + mu.0.a0w * (1 - gamma.a0)
  
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
  IF.0pw <- mu.a0 + 
    ifelse(dat.sim$A == 0, 
           ((dat.sim$R + dat.sim$S / eta.1) * dat.sim$Y - mu.a0 +
              (1 - dat.sim$R) * (1 - dat.sim$S / eta.1) * mu.0.a0) / (1 - pi.1w),
           0)
  IF.1pw <- mu.a1 + 
    ifelse(dat.sim$A == 1, 
           ((dat.sim$R + dat.sim$S / eta.1) * dat.sim$Y - mu.a1 +
              (1 - dat.sim$R) * (1 - dat.sim$S / eta.1) * mu.0.a1) / pi.1w,
           0)
  IF.0mw <- mu.a0w + 
    ifelse(dat.sim$A == 0, 
           ((dat.sim$R + dat.sim$S / eta.1) * dat.sim$Y - mu.a0w +
              (1 - dat.sim$R) * (1 - dat.sim$S / eta.1) * mu.0.a0w) / (1 - pi.1),
           0)
  IF.1mw <- mu.a1w + 
    ifelse(dat.sim$A == 1, 
           ((dat.sim$R + dat.sim$S / eta.1) * dat.sim$Y - mu.a1w +
              (1 - dat.sim$R) * (1 - dat.sim$S / eta.1) * mu.0.a1w) / pi.1,
           0)
  IF.0bw <- mu.a0w + 
    ifelse(dat.sim$A == 0, 
           ((dat.sim$R + dat.sim$S / eta.1) * dat.sim$Y - mu.a0w +
              (1 - dat.sim$R) * (1 - dat.sim$S / eta.1) * mu.0.a0w) / (1 - pi.1w),
           0)
  IF.1bw <- mu.a1w + 
    ifelse(dat.sim$A == 1, 
           ((dat.sim$R + dat.sim$S / eta.1) * dat.sim$Y - mu.a1w +
              (1 - dat.sim$R) * (1 - dat.sim$S / eta.1) * mu.0.a1w) / pi.1w,
           0)
  IPW.MAR.0 <- dat.sim$R * (1 - dat.sim$A) * dat.sim$Y /
    ((1 - pi.1) * gamma.a0)
  IPW.MAR.1 <- dat.sim$R * dat.sim$A * dat.sim$Y /
    (pi.1 * gamma.a1)
  IF.MAR.0 <- mu.1.a0 + dat.sim$R * (1 - dat.sim$A) * (dat.sim$Y - mu.1.a0) /
    ((1 - pi.1) * gamma.a0)
  IF.MAR.1 <- mu.1.a1 + dat.sim$R * dat.sim$A * (dat.sim$Y - mu.1.a1) /
    (pi.1 * gamma.a1)
  IF.MAR.eff.0 <- mu.e.a0 + (dat.sim$R + dat.sim$S) * (1 - dat.sim$A) *
    (dat.sim$Y - mu.e.a0) / ((1 - pi.1) * (gamma.a0 + (1 - gamma.a0) * eta.1))
  IF.MAR.eff.1 <- mu.e.a1 + (dat.sim$R + dat.sim$S) * dat.sim$A *
    (dat.sim$Y - mu.e.a1) / (pi.1 * (gamma.a1 + (1 - gamma.a1) * eta.1))
  
  ## for computing confidence intervals and Dominik Rothenhausler's criterion
  sigma.0.2[r] <- mean((IF.1 - IF.0 - mean(IF.1 - IF.0))^2)
  sigma.1.2[r] <- mean((IF.MAR.eff.1 - IF.MAR.eff.0 - 
                          mean(IF.MAR.eff.1 - IF.MAR.eff.0))^2)
  tau.1.2[r] <- mean((IF.1 - IF.0 - IF.MAR.eff.1 + IF.MAR.eff.0 - 
                        mean(IF.1 - IF.0 - IF.MAR.eff.1 + IF.MAR.eff.0))^2)
  
  res[r,] <- c(mean(IF.0), mean(IF.1),
               mean(IF.0pw), mean(IF.1pw),
               mean(IF.0mw), mean(IF.1mw),
               mean(IF.0bw), mean(IF.1bw),
               mean(mu.1.a0), mean(mu.1.a1),
               mean(IPW.MAR.0), mean(IPW.MAR.1),
               mean(IF.MAR.0), mean(IF.MAR.1),
               mean(IF.MAR.eff.0), mean(IF.MAR.eff.1))
}

res.long <- cbind.data.frame(c(res[,2] - res[,1],
                               res[,4] - res[,3],
                               res[,6] - res[,5],
                               res[,8] - res[,7],
                               res[,10] - res[,9],
                               res[,12] - res[,11],
                               res[,14] - res[,13],
                               res[,16] - res[,15]))
colnames(res.long) <- "value"
res.long$est <- c(rep("IF-DS", R),
                  rep("IF-DS-PW", R),
                  rep("IF-DS-MW", R),
                  rep("IF-DS-BW", R),
                  rep("OR-MAR", R),
                  rep("IPW-MAR", R),
                  rep("IF-MAR", R),
                  rep("IF-MAR-EFF", R))
## compute Dominik Rothenhausler criterion
R.hat <- cbind(sigma.0.2 / n, (sigma.1.2 - tau.1.2)/n + 
                 ((res[,16] - res[,15]) - (res[,2] - res[,1]))^2)
R.mod <- cbind(sigma.0.2 / n, sigma.1.2/n + 
                 pmax(((res[,16] - res[,15]) - (res[,2] - res[,1]))^2 - tau.1.2/n, 0))
res.long <- 
  rbind.data.frame(res.long,
                   cbind.data.frame(value = ifelse(R.hat[,2] < R.hat[,1],
                                                   res[,16] - res[,15],
                                                   res[,2] - res[,1]),
                                    est = rep("IF-adapt", R)),
                   cbind.data.frame(value = ifelse(R.mod[,2] < R.mod[,1],
                                                   res[,16] - res[,15],
                                                   res[,2] - res[,1]),
                                    est = rep("IF-adapt-mod", R)))
res.long$est <- factor(res.long$est,
                       levels = c("IF-DS","IF-DS-PW","IF-DS-MW","IF-DS-BW",
                                  "IF-MAR","OR-MAR","IPW-MAR","IF-MAR-EFF",
                                  "IF-adapt", "IF-adapt-mod"),
                       ordered = T)
# save(res.long, file = "RES_long_000.Rdata")
# load("RES_long_016.RData")
# res.long16 <- res.long

plot00 <- ggplot(res.long00[1:40000,], aes(x = factor(est), y = value)) +
  geom_boxplot() +
  geom_hline(yintercept = chi.1 - chi.0, linetype = "dashed", color = "red") +
  xlab("Estimator") +
  ylab("Estimated ATE")

plot16 <- ggplot(res.long16[1:40000,], aes(x = factor(est), y = value)) +
  geom_boxplot() +
  geom_hline(yintercept = chi.1 - chi.0, linetype = "dashed", color = "red") +
  xlab("Estimator") +
  ylab("Estimated ATE")

plot32 <- ggplot(res.long32[1:40000,], aes(x = factor(est), y = value)) +
  geom_boxplot() +
  geom_hline(yintercept = chi.1 - chi.0, linetype = "dashed", color = "red") +
  xlab("Estimator") +
  ylab("Estimated ATE")

library(ggpubr)
ggarrange(plot00, plot16, plot32, labels = c("a)","b)","c)"), ncol=1, nrow=3)

# theta_g = theta_0
xs <- seq((0 - 1) * (n * (var(res.long$value[res.long$est == "IF-DS"] - 
                                res.long$value[res.long$est == "IF-MAR-EFF"]))),
          3, by = 0.01)
ys <- dchisq(xs / (n * (var(res.long$value[res.long$est == "IF-DS"] -
                              res.long$value[res.long$est == "IF-MAR-EFF"]))) + 1, 1) / 
  (n * (var(res.long$value[res.long$est == "IF-DS"] - 
       res.long$value[res.long$est == "IF-MAR-EFF"])))
hist(n*(R.hat[,2] - var(res.long$value[res.long$est == "IF-MAR-EFF"])),
     breaks = 70, freq = F, main = expression('Distribution of' ~ hat(R)),
     xlab = expression(n(hat(R)(1) - MSE(1))))
lines(xs, ys, col = 'red')

# theta_g != theta_0
thresh <- max(abs(sqrt(n)*(R.hat[,2] - 
                             mean((res.long$value[res.long$est == "IF-MAR-EFF"] - 
                                     (chi.1 - chi.0))^2)))) * 1.05
xs <- seq(-thresh, thresh, by = thresh/1000)
ys <- dnorm(xs, 0, 2 * abs(mean(res.long$value[res.long$est == "IF-MAR-EFF"] - 
                                   (chi.1 - chi.0))) * 
              sqrt(n * (var(res.long$value[res.long$est == "IF-DS"] -
                              res.long$value[res.long$est == "IF-MAR-EFF"]))))
hist(sqrt(n)*(R.hat[,2] - mean((res.long$value[res.long$est == "IF-MAR-EFF"] - 
                                  (chi.1 - chi.0))^2)),
     breaks = 40, freq = F, main = expression('Distribution of' ~ hat(R)),
     xlab = expression(sqrt(n)(hat(R)(1) - MSE(1))))
lines(xs, ys, col = 'red')

## MSEs
mean((res.long$value[res.long$est == "IF-DS"] - (chi.1 - chi.0))^2)
mean((res.long$value[res.long$est == "IF-MAR-EFF"] - (chi.1 - chi.0))^2)
mean((res.long$value[res.long$est == "IF-adapt"] - (chi.1 - chi.0))^2)
mean((res.long$value[res.long$est == "IF-adapt-mod"] - (chi.1 - chi.0))^2)

## Variance estimators
hist(sigma.0.2 / sqrt(n))
hist(sigma.1.2 / sqrt(n))
hist(tau.1.2 / sqrt(n))


