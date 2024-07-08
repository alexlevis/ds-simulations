library(tidyverse)
library(SuperLearner)
library(mvtnorm)

# setwd("~/double_sampling/")
# dat <- read.csv("cases_yr3.csv")
# dat$R <- 1 * (! is.na(dat$bmi_change))
# mean(1 - dat$R) ## marginal missingness probability
# hist(dat$bmi_change, xlab = "BMI at 3 years - BMI at baseline", main = "Change in BMI")
# hist(dat$bmi_change[dat$bs_type == "SLEEVE"], xlab = "BMI at 3 years - BMI at baseline", main = "Change in BMI for Sleeve")
# hist(dat$bmi_change[dat$bs_type == "RYGB"], xlab = "BMI at 3 years - BMI at baseline", main = "Change in BMI for RYGB")
# dat <- dat[!is.na(dat$bmi_base),]      # need baseline BMI
# dat <- dat[!(dat$inscat == "Other"),]  # rule out rarest insurance type
# dat$mh_cat2 <- as.factor(ifelse(dat$mh_cat == "5: No preop mental illness", "None",
#                                 ifelse(dat$mh_cat == "3: Mild-to-Moderation Anxiety/Depression",
#                                        "Mild-Moderate Anxiety/Depression", "Other")))
# 
# ## MAR missingness mechanism
# ## R ~ A + X
# set.seed(498)
# mar.mod <- SuperLearner(dat$R, dat[,c("bs_type", "site", "year", "age", "gender","raceeth",
#                                       "util_count", "insulin", "elix_cat", "htn_dx",
#                                       "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
#                                       "diabetes", "days_ip", "dysl", "statins", "nonstatins",
#                                       "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
#                                       "retinopathy", "neuropathy", "bmi_base")],
#                         family = "binomial",
#                         SL.library = c("SL.glm", "SL.ranger","SL.rpart"))

# dat$MAR_probs <- pmax(pmin(mar.mod$SL.predict, 0.95), 0.05)
# write.csv(dat, file = "cases_yr3_MAR.csv", row.names = F)
dat <- read.csv("cases_yr3_MAR.csv")

## percent total weight change, as in Arterburn paper
dat$ptwc <- 100 * (dat$bmi_change / dat$bmi_base)
## standardize BMI ?
dat$bmi_change_std <- (dat$bmi_change - mean(dat$bmi_change, na.rm = T)) / 
  sd(dat$bmi_change, na.rm = T)
hist(dat$MAR_probs, xlab = "Missingness Probabilities", main = "")
dat$MNAR_probs <- pmax(pmin(plogis(qlogis(dat$MAR_probs) + 
                           0.7 * dat$bmi_change_std + 
                           0.8 * dat$bmi_change_std * (1 * (dat$bs_type == "SLEEVE")) +
                           -1.2 * dat$bmi_change_std * dat$diabetes +
                           0.6 * dat$bmi_change_std * (1 * (dat$inscat != "Commercial"))), 0.95), 0.05)

## three datasets
dat.cc <- dat[dat$R == 1, ]
dat.mar <- dat[dat$R == 1, ]
dat.mnar <- dat[dat$R == 1, ]

set.seed(514)
chunk <- 500

dat.mar$R <- rbinom(nrow(dat.mar), 1, dat.mar$MAR_probs)
S.mar.remain <- (1:nrow(dat.mar))[dat.mar$R == 0]
S.mar.0500 <- sample(S.mar.remain, chunk, replace = F); S.mar.remain <- S.mar.remain[! S.mar.remain %in% S.mar.0500]
S.mar.1000 <- sample(S.mar.remain, chunk, replace = F); S.mar.remain <- S.mar.remain[! S.mar.remain %in% S.mar.1000]
S.mar.1500 <- sample(S.mar.remain, chunk, replace = F)
dat.mar$S1 <- 1 * ((1:nrow(dat.mar)) %in% S.mar.0500)
dat.mar$S2 <- 1 * ((1:nrow(dat.mar)) %in% c(S.mar.0500, S.mar.1000))
dat.mar$S3 <- 1 * ((1:nrow(dat.mar)) %in% c(S.mar.0500, S.mar.1000, S.mar.1500))
dat.mar$ptwc3 <- dat.mar$ptwc2 <- dat.mar$ptwc1 <- dat.mar$ptwc
dat.mar$ptwc1[dat.mar$R == 0 & dat.mar$S1 == 0] <- NA
dat.mar$ptwc2[dat.mar$R == 0 & dat.mar$S2 == 0] <- NA
dat.mar$ptwc3[dat.mar$R == 0 & dat.mar$S3 == 0] <- NA
eta.0.mar.1 <- sum(dat.mar$S1) / sum(dat.mar$R == 0)
eta.0.mar.2 <- sum(dat.mar$S2) / sum(dat.mar$R == 0)
eta.0.mar.3 <- sum(dat.mar$S3) / sum(dat.mar$R == 0)

dat.mnar$R <- rbinom(nrow(dat.mnar), 1, dat.mnar$MNAR_probs)
S.mnar.remain <- (1:nrow(dat.mnar))[dat.mnar$R == 0]
S.mnar.0500 <- sample(S.mnar.remain, chunk, replace = F); S.mnar.remain <- S.mnar.remain[! S.mnar.remain %in% S.mnar.0500]
S.mnar.1000 <- sample(S.mnar.remain, chunk, replace = F); S.mnar.remain <- S.mnar.remain[! S.mnar.remain %in% S.mnar.1000]
S.mnar.1500 <- sample(S.mnar.remain, chunk, replace = F)
dat.mnar$S1 <- 1 * ((1:nrow(dat.mnar)) %in% S.mnar.0500)
dat.mnar$S2 <- 1 * ((1:nrow(dat.mnar)) %in% c(S.mnar.0500, S.mnar.1000))
dat.mnar$S3 <- 1 * ((1:nrow(dat.mnar)) %in% c(S.mnar.0500, S.mnar.1000, S.mnar.1500))
dat.mnar$ptwc3 <- dat.mnar$ptwc2 <- dat.mnar$ptwc1 <- dat.mnar$ptwc
dat.mnar$ptwc1[dat.mnar$R == 0 & dat.mnar$S1 == 0] <- NA
dat.mnar$ptwc2[dat.mnar$R == 0 & dat.mnar$S2 == 0] <- NA
dat.mnar$ptwc3[dat.mnar$R == 0 & dat.mnar$S3 == 0] <- NA
eta.0.mnar.1 <- sum(dat.mnar$S1) / sum(dat.mnar$R == 0)
eta.0.mnar.2 <- sum(dat.mnar$S2) / sum(dat.mnar$R == 0)
eta.0.mnar.3 <- sum(dat.mnar$S3) / sum(dat.mnar$R == 0)

### Analysis with full data
dat.cc$A <- 1 * (dat.cc$bs_type == "RYGB") # ATE will be RYGB - SLEEVE
out.cc <- SuperLearner(dat.cc$ptwc, 
                       dat.cc[,c("A", "site", "year", "age", "gender","raceeth",
                                 "util_count", "insulin", "elix_cat", "htn_dx",
                                 "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                 "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                 "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                 "retinopathy", "neuropathy", "bmi_base")],
                       SL.library = c("SL.glm", "SL.ranger","SL.rpart"))
prop.cc <- SuperLearner(dat.cc$A, dat.cc[,c("site", "year", "age", "gender","raceeth",
                                            "util_count", "insulin", "elix_cat", "htn_dx",
                                            "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                            "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                            "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                            "retinopathy", "neuropathy", "bmi_base")],
                        family = "binomial",
                        SL.library = c("SL.glm", "SL.ranger","SL.rpart"))

## nuisance function estimates
dat.cc.1 <- cbind.data.frame(dat.cc[,-42], A = 1)
dat.cc.0 <- cbind.data.frame(dat.cc[,-42], A = 0)
mu.cc.1 <- predict.SuperLearner(out.cc, dat.cc.1)$pred
mu.cc.0 <- predict.SuperLearner(out.cc, dat.cc.0)$pred

pi.cc <- pmin(0.98, pmax(0.02, predict.SuperLearner(prop.cc, type = "response")$pred))
# pi.cc <- predict(prop.cc, type = "response")

## IF and final estimator
IF.cc <- mu.cc.1 - mu.cc.0 + dat.cc$A * (dat.cc$ptwc - mu.cc.1) / pi.cc - 
  (1 - dat.cc$A) * (dat.cc$ptwc - mu.cc.0) / (1 - pi.cc)

chi.cc <- mean(IF.cc)
chi.var.cc <- mean((IF.cc - chi.cc)^2)
chi.cc                                                           # point estimate
chi.cc + c(-1, 1)*qnorm(0.975) * sqrt(chi.var.cc / nrow(dat.cc)) # 95% CI


### MAR analyses
dat.mar$A <- 1 * (dat.mar$bs_type == "RYGB") # ATE will be RYGB - SLEEVE
out.mar <- SuperLearner(dat.mar$ptwc[dat.mar$R == 1], 
                       dat.mar[dat.mar$R == 1,
                              c("A", "site", "year", "age", "gender","raceeth",
                                 "util_count", "insulin", "elix_cat", "htn_dx",
                                 "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                 "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                 "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                 "retinopathy", "neuropathy", "bmi_base")],
                       SL.library = c("SL.glm", "SL.ranger","SL.rpart"))
miss.mar <- SuperLearner(dat.mar$R, dat.mar[,c("A", "site", "year", "age", "gender","raceeth",
                                               "util_count", "insulin", "elix_cat", "htn_dx",
                                               "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                               "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                               "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                               "retinopathy", "neuropathy", "bmi_base")],
                         family = "binomial",
                         SL.library = c("SL.glm", "SL.ranger","SL.rpart"))

out.mar.S1 <- SuperLearner(dat.mar$ptwc[dat.mar$S1 == 1],
                           dat.mar[dat.mar$S1 == 1,
                                   c("A", "site", "year", "age", "gender","raceeth",
                                    "util_count", "insulin", "elix_cat", "htn_dx",
                                    "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                    "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                    "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                    "retinopathy", "neuropathy", "bmi_base")],
                          SL.library = c("SL.glm", "SL.ranger", "SL.rpart"))
out.mar.comb1 <- SuperLearner(dat.mar$ptwc[dat.mar$R == 1 | dat.mar$S1 == 1], 
                             dat.mar[dat.mar$R == 1 | dat.mar$S1 == 1,
                                     c("A", "site", "year", "age", "gender","raceeth",
                                       "util_count", "insulin", "elix_cat", "htn_dx",
                                       "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                       "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                       "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                       "retinopathy", "neuropathy", "bmi_base")],
                             SL.library = c("SL.glm", "SL.ranger","SL.rpart"))

out.mar.S2 <- SuperLearner(dat.mar$ptwc[dat.mar$S2 == 1],
                           dat.mar[dat.mar$S2 == 1,
                                   c("A", "site", "year", "age", "gender","raceeth",
                                     "util_count", "insulin", "elix_cat", "htn_dx",
                                     "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                     "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                     "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                     "retinopathy", "neuropathy", "bmi_base")],
                           SL.library = c("SL.glm", "SL.ranger", "SL.rpart"))
out.mar.comb2 <- SuperLearner(dat.mar$ptwc[dat.mar$R == 1 | dat.mar$S2 == 1], 
                              dat.mar[dat.mar$R == 1 | dat.mar$S2 == 1,
                                      c("A", "site", "year", "age", "gender","raceeth",
                                        "util_count", "insulin", "elix_cat", "htn_dx",
                                        "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                        "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                        "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                        "retinopathy", "neuropathy", "bmi_base")],
                              SL.library = c("SL.glm", "SL.ranger","SL.rpart"))

out.mar.S3 <- SuperLearner(dat.mar$ptwc[dat.mar$S3 == 1],
                           dat.mar[dat.mar$S3 == 1,
                                   c("A", "site", "year", "age", "gender","raceeth",
                                     "util_count", "insulin", "elix_cat", "htn_dx",
                                     "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                     "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                     "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                     "retinopathy", "neuropathy", "bmi_base")],
                           SL.library = c("SL.glm", "SL.ranger", "SL.rpart"))
out.mar.comb3 <- SuperLearner(dat.mar$ptwc[dat.mar$R == 1 | dat.mar$S3 == 1], 
                              dat.mar[dat.mar$R == 1 | dat.mar$S3 == 1,
                                      c("A", "site", "year", "age", "gender","raceeth",
                                        "util_count", "insulin", "elix_cat", "htn_dx",
                                        "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                        "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                        "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                        "retinopathy", "neuropathy", "bmi_base")],
                              SL.library = c("SL.glm", "SL.ranger","SL.rpart"))

## nuisance function estimates
dat.mar.1 <- cbind.data.frame(dat.mar[,-48], A = 1)
dat.mar.0 <- cbind.data.frame(dat.mar[,-48], A = 0)
mu.mar.1 <- predict.SuperLearner(out.mar, dat.mar.1)$pred
mu.mar.0 <- predict.SuperLearner(out.mar, dat.mar.0)$pred
gamma.mar <- predict.SuperLearner(miss.mar, dat.mar, type = "response")$pred
gamma.mar.1 <- predict.SuperLearner(miss.mar, dat.mar.1, type = "response")$pred
gamma.mar.0 <- predict.SuperLearner(miss.mar, dat.mar.0, type = "response")$pred

mu.mar.1.S1 <- predict.SuperLearner(out.mar.S1, dat.mar.1)$pred
mu.mar.0.S1 <- predict.SuperLearner(out.mar.S1, dat.mar.0)$pred
mu.mar.1.S2 <- predict.SuperLearner(out.mar.S2, dat.mar.1)$pred
mu.mar.0.S2 <- predict.SuperLearner(out.mar.S2, dat.mar.0)$pred
mu.mar.1.S3 <- predict.SuperLearner(out.mar.S3, dat.mar.1)$pred
mu.mar.0.S3 <- predict.SuperLearner(out.mar.S3, dat.mar.0)$pred

mu.mar.eff1.1 <- predict.SuperLearner(out.mar.comb1, dat.mar.1)$pred
mu.mar.eff1.0 <- predict.SuperLearner(out.mar.comb1, dat.mar.0)$pred
mu.mar.eff2.1 <- predict.SuperLearner(out.mar.comb2, dat.mar.1)$pred
mu.mar.eff2.0 <- predict.SuperLearner(out.mar.comb2, dat.mar.0)$pred
mu.mar.eff3.1 <- predict.SuperLearner(out.mar.comb3, dat.mar.1)$pred
mu.mar.eff3.0 <- predict.SuperLearner(out.mar.comb3, dat.mar.0)$pred

## IF and MAR-based incomplete data estimator
IF.mar <- mu.mar.1 - mu.mar.0 + ifelse(dat.mar$R == 1, 
  (dat.mar$A * (dat.mar$ptwc - mu.mar.1) / pi.cc - 
  (1 - dat.mar$A) * (dat.mar$ptwc - mu.mar.0) / (1 - pi.cc)), 0) / gamma.mar

chi.mar <- mean(IF.mar)
chi.var.mar <- mean((IF.mar - chi.mar)^2)
chi.mar                                                             # point estimate
chi.mar + c(-1, 1)*qnorm(0.975) * sqrt(chi.var.mar / nrow(dat.mar)) # 95% CI

## IF and nonparametric efficient estimator
mu.mar.np1.1 <- gamma.mar.1 * mu.mar.1 + (1 - gamma.mar.1) * mu.mar.1.S1
mu.mar.np1.0 <- gamma.mar.0 * mu.mar.0 + (1 - gamma.mar.0) * mu.mar.0.S1
IF.mar.np1 <- mu.mar.np1.1 - mu.mar.np1.0 + 
  dat.mar$A * ((dat.mar$R + dat.mar$S1 / eta.0.mar.1) * ifelse(dat.mar$R == 1 | dat.mar$S1 == 1,
                                                        dat.mar$ptwc, 0) - 
                 mu.mar.np1.1 +
                 (1 - dat.mar$R) * (1 - dat.mar$S1 / eta.0.mar.1) * mu.mar.1.S1) / pi.cc - 
  (1 - dat.mar$A) * ((dat.mar$R + dat.mar$S1 / eta.0.mar.1) * ifelse(dat.mar$R == 1 | dat.mar$S1 == 1, 
                                                              dat.mar$ptwc, 0) - 
                       mu.mar.np1.0 +
                       (1 - dat.mar$R) * (1 - dat.mar$S1 / eta.0.mar.1) * mu.mar.0.S1) / (1 - pi.cc)

tau.mar1 <- mean(IF.mar.np1)
tau.var.mar1 <- mean((IF.mar.np1 - tau.mar1)^2)
tau.mar1                                                             # point estimate
tau.mar1 + c(-1, 1)*qnorm(0.975) * sqrt(tau.var.mar1 / nrow(dat.mar)) # 95% CI

mu.mar.np2.1 <- gamma.mar.1 * mu.mar.1 + (1 - gamma.mar.1) * mu.mar.1.S2
mu.mar.np2.0 <- gamma.mar.0 * mu.mar.0 + (1 - gamma.mar.0) * mu.mar.0.S2
IF.mar.np2 <- mu.mar.np2.1 - mu.mar.np2.0 + 
  dat.mar$A * ((dat.mar$R + dat.mar$S2 / eta.0.mar.2) * ifelse(dat.mar$R == 1 | dat.mar$S2 == 1,
                                                               dat.mar$ptwc, 0) - 
                 mu.mar.np2.1 +
                 (1 - dat.mar$R) * (1 - dat.mar$S2 / eta.0.mar.2) * mu.mar.1.S2) / pi.cc - 
  (1 - dat.mar$A) * ((dat.mar$R + dat.mar$S2 / eta.0.mar.2) * ifelse(dat.mar$R == 1 | dat.mar$S2 == 1, 
                                                                     dat.mar$ptwc, 0) - 
                       mu.mar.np2.0 +
                       (1 - dat.mar$R) * (1 - dat.mar$S2 / eta.0.mar.2) * mu.mar.0.S2) / (1 - pi.cc)

tau.mar2 <- mean(IF.mar.np2)
tau.var.mar2 <- mean((IF.mar.np2 - tau.mar2)^2)
tau.mar2                                                             # point estimate
tau.mar2 + c(-1, 1)*qnorm(0.975) * sqrt(tau.var.mar2 / nrow(dat.mar)) # 95% CI

mu.mar.np3.1 <- gamma.mar.1 * mu.mar.1 + (1 - gamma.mar.1) * mu.mar.1.S3
mu.mar.np3.0 <- gamma.mar.0 * mu.mar.0 + (1 - gamma.mar.0) * mu.mar.0.S3
IF.mar.np3 <- mu.mar.np3.1 - mu.mar.np3.0 + 
  dat.mar$A * ((dat.mar$R + dat.mar$S3 / eta.0.mar.3) * ifelse(dat.mar$R == 1 | dat.mar$S3 == 1,
                                                               dat.mar$ptwc, 0) - 
                 mu.mar.np3.1 +
                 (1 - dat.mar$R) * (1 - dat.mar$S3 / eta.0.mar.3) * mu.mar.1.S3) / pi.cc - 
  (1 - dat.mar$A) * ((dat.mar$R + dat.mar$S3 / eta.0.mar.3) * ifelse(dat.mar$R == 1 | dat.mar$S3 == 1, 
                                                                     dat.mar$ptwc, 0) - 
                       mu.mar.np3.0 +
                       (1 - dat.mar$R) * (1 - dat.mar$S3 / eta.0.mar.3) * mu.mar.0.S3) / (1 - pi.cc)

tau.mar3 <- mean(IF.mar.np3)
tau.var.mar3 <- mean((IF.mar.np3 - tau.mar3)^2)
tau.mar3                                                             # point estimate
tau.mar3 + c(-1, 1)*qnorm(0.975) * sqrt(tau.var.mar3 / nrow(dat.mar)) # 95% CI


## IF and MAR-based efficient estimator
IF.mar.eff1 <- mu.mar.eff1.1 - mu.mar.eff1.0 + 
  ifelse(dat.mar$R == 1 | dat.mar$S1 == 1, 
         (dat.mar$A * (dat.mar$ptwc - mu.mar.eff1.1) / pi.cc - 
            (1 - dat.mar$A) * (dat.mar$ptwc - mu.mar.eff1.0) / (1 - pi.cc)), 0) / 
  (1 - (1 - gamma.mar) * (1 - eta.0.mar.1))

taustar.mar1 <- mean(IF.mar.eff1)
taustar.var.mar1 <- mean((IF.mar.eff1 - taustar.mar1)^2)
taustar.mar1                                                                 # point estimate
taustar.mar1 + c(-1, 1)*qnorm(0.975) * sqrt(taustar.var.mar1 / nrow(dat.mar)) # 95% CI

IF.mar.eff2 <- mu.mar.eff2.1 - mu.mar.eff2.0 + 
  ifelse(dat.mar$R == 1 | dat.mar$S2 == 1, 
         (dat.mar$A * (dat.mar$ptwc - mu.mar.eff2.1) / pi.cc - 
            (1 - dat.mar$A) * (dat.mar$ptwc - mu.mar.eff2.0) / (1 - pi.cc)), 0) / 
  (1 - (1 - gamma.mar) * (1 - eta.0.mar.2))

taustar.mar2 <- mean(IF.mar.eff2)
taustar.var.mar2 <- mean((IF.mar.eff2 - taustar.mar2)^2)
taustar.mar2                                                                 # point estimate
taustar.mar2 + c(-1, 1)*qnorm(0.975) * sqrt(taustar.var.mar2 / nrow(dat.mar)) # 95% CI

IF.mar.eff3 <- mu.mar.eff3.1 - mu.mar.eff3.0 + 
  ifelse(dat.mar$R == 1 | dat.mar$S3 == 1, 
         (dat.mar$A * (dat.mar$ptwc - mu.mar.eff3.1) / pi.cc - 
            (1 - dat.mar$A) * (dat.mar$ptwc - mu.mar.eff3.0) / (1 - pi.cc)), 0) / 
  (1 - (1 - gamma.mar) * (1 - eta.0.mar.3))

taustar.mar3 <- mean(IF.mar.eff3)
taustar.var.mar3 <- mean((IF.mar.eff3 - taustar.mar3)^2)
taustar.mar3                                                                 # point estimate
taustar.mar3 + c(-1, 1)*qnorm(0.975) * sqrt(taustar.var.mar3 / nrow(dat.mar)) # 95% CI


## adaptive estimator
diff.var.mar1 <- mean(((IF.mar.np1 - IF.mar.eff1) - (tau.mar1 - taustar.mar1))^2)
R.mod.mar1 <- c(tau.var.mar1 / nrow(dat.mar), taustar.var.mar1 / nrow(dat.mar) + 
             max((tau.mar1 - taustar.mar1)^2 - diff.var.mar1 / nrow(dat.mar), 0))
cross.cov.mar1 <- mean((IF.mar.np1 - tau.mar1) * (IF.mar.eff1 -  taustar.mar1))
tauadapt.mar1 <- ifelse(R.mod.mar1[1] < R.mod.mar1[2], tau.mar1, taustar.mar1)

diff.var.mar2 <- mean(((IF.mar.np2 - IF.mar.eff2) - (tau.mar2 - taustar.mar2))^2)
R.mod.mar2 <- c(tau.var.mar2 / nrow(dat.mar), taustar.var.mar2 / nrow(dat.mar) + 
                  max((tau.mar2 - taustar.mar2)^2 - diff.var.mar2 / nrow(dat.mar), 0))
cross.cov.mar2 <- mean((IF.mar.np2 - tau.mar2) * (IF.mar.eff2 -  taustar.mar2))
tauadapt.mar2 <- ifelse(R.mod.mar2[1] < R.mod.mar2[2], tau.mar2, taustar.mar2)

diff.var.mar3 <- mean(((IF.mar.np3 - IF.mar.eff3) - (tau.mar3 - taustar.mar3))^2)
R.mod.mar3 <- c(tau.var.mar3 / nrow(dat.mar), taustar.var.mar3 / nrow(dat.mar) + 
                  max((tau.mar3 - taustar.mar3)^2 - diff.var.mar3 / nrow(dat.mar), 0))
cross.cov.mar3 <- mean((IF.mar.np3 - tau.mar3) * (IF.mar.eff3 -  taustar.mar3))
tauadapt.mar3 <- ifelse(R.mod.mar3[1] < R.mod.mar3[2], tau.mar3, taustar.mar3)

## compute CIs for adaptive estimator
set.seed(264)
M <- 5000
beta.grid <- seq(-300, 300, by = 0.01)

Sigma1 <- rbind(c(tau.var.mar1, cross.cov.mar1),
               c(cross.cov.mar1, taustar.var.mar1))
Z <- rmvnorm(M, sqrt(nrow(dat.mar)) * c(tau.mar1, taustar.mar1), Sigma1)
## compute indicator of the event A_0^{est}
A.g <- 1 * (tau.var.mar1 < pmax((Z[,2] - Z[,1])^2, 0) + taustar.var.mar1)
b.beta <- sapply(beta.grid, function(x) {
  mean((Z[,1] - sqrt(nrow(dat.mar)) * tau.mar1 <= x) & A.g == 1) +
    mean((Z[,2] - sqrt(nrow(dat.mar)) * tau.mar1 <= x) & A.g == 0)
}, simplify = 0)
b.025 <- beta.grid[which.min(abs(b.beta - 0.025))]
b.975 <- beta.grid[which.min(abs(b.beta - 0.975))]
tauadapt.CI.mar1 <- c(tauadapt.mar1 - b.975 / sqrt(nrow(dat.mar)),
                     tauadapt.mar1 - b.025 / sqrt(nrow(dat.mar)))

Sigma2 <- rbind(c(tau.var.mar2, cross.cov.mar2),
                c(cross.cov.mar2, taustar.var.mar2))
Z <- rmvnorm(M, sqrt(nrow(dat.mar)) * c(tau.mar2, taustar.mar2), Sigma2)
## compute indicator of the event A_0^{est}
A.g <- 1 * (tau.var.mar2 < pmax((Z[,2] - Z[,1])^2, 0) + taustar.var.mar2)
b.beta <- sapply(beta.grid, function(x) {
  mean((Z[,1] - sqrt(nrow(dat.mar)) * tau.mar2 <= x) & A.g == 1) +
    mean((Z[,2] - sqrt(nrow(dat.mar)) * tau.mar2 <= x) & A.g == 0)
}, simplify = 0)
b.025 <- beta.grid[which.min(abs(b.beta - 0.025))]
b.975 <- beta.grid[which.min(abs(b.beta - 0.975))]
tauadapt.CI.mar2 <- c(tauadapt.mar2 - b.975 / sqrt(nrow(dat.mar)),
                      tauadapt.mar2 - b.025 / sqrt(nrow(dat.mar)))

Sigma3 <- rbind(c(tau.var.mar3, cross.cov.mar3),
                c(cross.cov.mar3, taustar.var.mar3))
Z <- rmvnorm(M, sqrt(nrow(dat.mar)) * c(tau.mar3, taustar.mar3), Sigma3)
## compute indicator of the event A_0^{est}
A.g <- 1 * (tau.var.mar3 < pmax((Z[,2] - Z[,1])^2, 0) + taustar.var.mar3)
b.beta <- sapply(beta.grid, function(x) {
  mean((Z[,1] - sqrt(nrow(dat.mar)) * tau.mar3 <= x) & A.g == 1) +
    mean((Z[,2] - sqrt(nrow(dat.mar)) * tau.mar3 <= x) & A.g == 0)
}, simplify = 0)
b.025 <- beta.grid[which.min(abs(b.beta - 0.025))]
b.975 <- beta.grid[which.min(abs(b.beta - 0.975))]
tauadapt.CI.mar3 <- c(tauadapt.mar3 - b.975 / sqrt(nrow(dat.mar)),
                      tauadapt.mar3 - b.025 / sqrt(nrow(dat.mar)))

### MNAR analyses
dat.mnar$A <- 1 * (dat.mnar$bs_type == "RYGB") # ATE will be RYGB - SLEEVE
out.mnar <- SuperLearner(dat.mnar$ptwc[dat.mnar$R == 1], 
                        dat.mnar[dat.mnar$R == 1,
                                c("A", "site", "year", "age", "gender","raceeth",
                                  "util_count", "insulin", "elix_cat", "htn_dx",
                                  "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                  "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                  "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                  "retinopathy", "neuropathy", "bmi_base")],
                        SL.library = c("SL.glm", "SL.ranger","SL.rpart"))
miss.mnar <- SuperLearner(dat.mnar$R, dat.mnar[,c("A", "site", "year", "age", "gender","raceeth",
                                               "util_count", "insulin", "elix_cat", "htn_dx",
                                               "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                               "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                               "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                               "retinopathy", "neuropathy", "bmi_base")],
                         family = "binomial",
                         SL.library = c("SL.glm", "SL.ranger","SL.rpart"))

out.mnar.S1 <- SuperLearner(dat.mnar$ptwc[dat.mnar$S1 == 1],
                           dat.mnar[dat.mnar$S1 == 1,
                                   c("A", "site", "year", "age", "gender","raceeth",
                                     "util_count", "insulin", "elix_cat", "htn_dx",
                                     "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                     "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                     "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                     "retinopathy", "neuropathy", "bmi_base")],
                           SL.library = c("SL.glm", "SL.ranger", "SL.rpart"))
out.mnar.comb1 <- SuperLearner(dat.mnar$ptwc[dat.mnar$R == 1 | dat.mnar$S1 == 1], 
                              dat.mnar[dat.mnar$R == 1 | dat.mnar$S1 == 1,
                                      c("A", "site", "year", "age", "gender","raceeth",
                                        "util_count", "insulin", "elix_cat", "htn_dx",
                                        "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                        "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                        "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                        "retinopathy", "neuropathy", "bmi_base")],
                              SL.library = c("SL.glm", "SL.ranger","SL.rpart"))

out.mnar.S2 <- SuperLearner(dat.mnar$ptwc[dat.mnar$S2 == 1],
                           dat.mnar[dat.mnar$S2 == 1,
                                   c("A", "site", "year", "age", "gender","raceeth",
                                     "util_count", "insulin", "elix_cat", "htn_dx",
                                     "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                     "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                     "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                     "retinopathy", "neuropathy", "bmi_base")],
                           SL.library = c("SL.glm", "SL.ranger", "SL.rpart"))
out.mnar.comb2 <- SuperLearner(dat.mnar$ptwc[dat.mnar$R == 1 | dat.mnar$S2 == 1], 
                              dat.mnar[dat.mnar$R == 1 | dat.mnar$S2 == 1,
                                      c("A", "site", "year", "age", "gender","raceeth",
                                        "util_count", "insulin", "elix_cat", "htn_dx",
                                        "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                        "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                        "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                        "retinopathy", "neuropathy", "bmi_base")],
                              SL.library = c("SL.glm", "SL.ranger","SL.rpart"))

out.mnar.S3 <- SuperLearner(dat.mnar$ptwc[dat.mnar$S3 == 1],
                           dat.mnar[dat.mnar$S3 == 1,
                                   c("A", "site", "year", "age", "gender","raceeth",
                                     "util_count", "insulin", "elix_cat", "htn_dx",
                                     "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                     "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                     "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                     "retinopathy", "neuropathy", "bmi_base")],
                           SL.library = c("SL.glm", "SL.ranger", "SL.rpart"))
out.mnar.comb3 <- SuperLearner(dat.mnar$ptwc[dat.mnar$R == 1 | dat.mnar$S3 == 1], 
                              dat.mnar[dat.mnar$R == 1 | dat.mnar$S3 == 1,
                                      c("A", "site", "year", "age", "gender","raceeth",
                                        "util_count", "insulin", "elix_cat", "htn_dx",
                                        "htn_ace", "htn_arb", "htn_non_acearb", "inscat",
                                        "diabetes", "days_ip", "dysl", "statins", "nonstatins",
                                        "smoking2", "cerebro", "cardio", "cad", "mh_cat2",
                                        "retinopathy", "neuropathy", "bmi_base")],
                              SL.library = c("SL.glm", "SL.ranger","SL.rpart"))

## nuisance function estimates
dat.mnar.1 <- cbind.data.frame(dat.mnar[,-48], A = 1)
dat.mnar.0 <- cbind.data.frame(dat.mnar[,-48], A = 0)
mu.mnar.1 <- predict.SuperLearner(out.mnar, dat.mnar.1)$pred
mu.mnar.0 <- predict.SuperLearner(out.mnar, dat.mnar.0)$pred
gamma.mnar <- predict.SuperLearner(miss.mnar, dat.mnar, type = "response")$pred
gamma.mnar.1 <- predict.SuperLearner(miss.mnar, dat.mnar.1, type = "response")$pred
gamma.mnar.0 <- predict.SuperLearner(miss.mnar, dat.mnar.0, type = "response")$pred

mu.mnar.1.S1 <- predict.SuperLearner(out.mnar.S1, dat.mnar.1)$pred
mu.mnar.0.S1 <- predict.SuperLearner(out.mnar.S1, dat.mnar.0)$pred
mu.mnar.1.S2 <- predict.SuperLearner(out.mnar.S2, dat.mnar.1)$pred
mu.mnar.0.S2 <- predict.SuperLearner(out.mnar.S2, dat.mnar.0)$pred
mu.mnar.1.S3 <- predict.SuperLearner(out.mnar.S3, dat.mnar.1)$pred
mu.mnar.0.S3 <- predict.SuperLearner(out.mnar.S3, dat.mnar.0)$pred

mu.mnar.eff1.1 <- predict.SuperLearner(out.mnar.comb1, dat.mnar.1)$pred
mu.mnar.eff1.0 <- predict.SuperLearner(out.mnar.comb1, dat.mnar.0)$pred
mu.mnar.eff2.1 <- predict.SuperLearner(out.mnar.comb2, dat.mnar.1)$pred
mu.mnar.eff2.0 <- predict.SuperLearner(out.mnar.comb2, dat.mnar.0)$pred
mu.mnar.eff3.1 <- predict.SuperLearner(out.mnar.comb3, dat.mnar.1)$pred
mu.mnar.eff3.0 <- predict.SuperLearner(out.mnar.comb3, dat.mnar.0)$pred

## IF and MAR-based incomplete data estimator
IF.mnar <- mu.mnar.1 - mu.mnar.0 + ifelse(dat.mnar$R == 1, 
                                       (dat.mnar$A * (dat.mnar$ptwc - mu.mnar.1) / pi.cc - 
                                          (1 - dat.mnar$A) * (dat.mnar$ptwc - mu.mnar.0) / (1 - pi.cc)), 0) / gamma.mnar

chi.mnar <- mean(IF.mnar)
chi.var.mnar <- mean((IF.mnar - chi.mnar)^2)
chi.mnar                                                             # point estimate
chi.mnar + c(-1, 1)*qnorm(0.975) * sqrt(chi.var.mnar / nrow(dat.mnar)) # 95% CI

## IF and nonparametric efficient estimator
mu.mnar.np1.1 <- gamma.mnar.1 * mu.mnar.1 + (1 - gamma.mnar.1) * mu.mnar.1.S1
mu.mnar.np1.0 <- gamma.mnar.0 * mu.mnar.0 + (1 - gamma.mnar.0) * mu.mnar.0.S1
IF.mnar.np1 <- mu.mnar.np1.1 - mu.mnar.np1.0 + 
  dat.mnar$A * ((dat.mnar$R + dat.mnar$S1 / eta.0.mnar.1) * ifelse(dat.mnar$R == 1 | dat.mnar$S1 == 1,
                                                               dat.mnar$ptwc, 0) - 
                 mu.mnar.np1.1 +
                 (1 - dat.mnar$R) * (1 - dat.mnar$S1 / eta.0.mnar.1) * mu.mnar.1.S1) / pi.cc - 
  (1 - dat.mnar$A) * ((dat.mnar$R + dat.mnar$S1 / eta.0.mnar.1) * ifelse(dat.mnar$R == 1 | dat.mnar$S1 == 1, 
                                                                     dat.mnar$ptwc, 0) - 
                       mu.mnar.np1.0 +
                       (1 - dat.mnar$R) * (1 - dat.mnar$S1 / eta.0.mnar.1) * mu.mnar.0.S1) / (1 - pi.cc)

tau.mnar1 <- mean(IF.mnar.np1)
tau.var.mnar1 <- mean((IF.mnar.np1 - tau.mnar1)^2)
tau.mnar1                                                             # point estimate
tau.mnar1 + c(-1, 1)*qnorm(0.975) * sqrt(tau.var.mnar1 / nrow(dat.mnar)) # 95% CI

mu.mnar.np2.1 <- gamma.mnar.1 * mu.mnar.1 + (1 - gamma.mnar.1) * mu.mnar.1.S2
mu.mnar.np2.0 <- gamma.mnar.0 * mu.mnar.0 + (1 - gamma.mnar.0) * mu.mnar.0.S2
IF.mnar.np2 <- mu.mnar.np2.1 - mu.mnar.np2.0 + 
  dat.mnar$A * ((dat.mnar$R + dat.mnar$S2 / eta.0.mnar.2) * ifelse(dat.mnar$R == 1 | dat.mnar$S2 == 1,
                                                               dat.mnar$ptwc, 0) - 
                 mu.mnar.np2.1 +
                 (1 - dat.mnar$R) * (1 - dat.mnar$S2 / eta.0.mnar.2) * mu.mnar.1.S2) / pi.cc - 
  (1 - dat.mnar$A) * ((dat.mnar$R + dat.mnar$S2 / eta.0.mnar.2) * ifelse(dat.mnar$R == 1 | dat.mnar$S2 == 1, 
                                                                     dat.mnar$ptwc, 0) - 
                       mu.mnar.np2.0 +
                       (1 - dat.mnar$R) * (1 - dat.mnar$S2 / eta.0.mnar.2) * mu.mnar.0.S2) / (1 - pi.cc)

tau.mnar2 <- mean(IF.mnar.np2)
tau.var.mnar2 <- mean((IF.mnar.np2 - tau.mnar2)^2)
tau.mnar2                                                             # point estimate
tau.mnar2 + c(-1, 1)*qnorm(0.975) * sqrt(tau.var.mnar2 / nrow(dat.mnar)) # 95% CI

mu.mnar.np3.1 <- gamma.mnar.1 * mu.mnar.1 + (1 - gamma.mnar.1) * mu.mnar.1.S3
mu.mnar.np3.0 <- gamma.mnar.0 * mu.mnar.0 + (1 - gamma.mnar.0) * mu.mnar.0.S3
IF.mnar.np3 <- mu.mnar.np3.1 - mu.mnar.np3.0 + 
  dat.mnar$A * ((dat.mnar$R + dat.mnar$S3 / eta.0.mnar.3) * ifelse(dat.mnar$R == 1 | dat.mnar$S3 == 1,
                                                               dat.mnar$ptwc, 0) - 
                 mu.mnar.np3.1 +
                 (1 - dat.mnar$R) * (1 - dat.mnar$S3 / eta.0.mnar.3) * mu.mnar.1.S3) / pi.cc - 
  (1 - dat.mnar$A) * ((dat.mnar$R + dat.mnar$S3 / eta.0.mnar.3) * ifelse(dat.mnar$R == 1 | dat.mnar$S3 == 1, 
                                                                     dat.mnar$ptwc, 0) - 
                       mu.mnar.np3.0 +
                       (1 - dat.mnar$R) * (1 - dat.mnar$S3 / eta.0.mnar.3) * mu.mnar.0.S3) / (1 - pi.cc)

tau.mnar3 <- mean(IF.mnar.np3)
tau.var.mnar3 <- mean((IF.mnar.np3 - tau.mnar3)^2)
tau.mnar3                                                             # point estimate
tau.mnar3 + c(-1, 1)*qnorm(0.975) * sqrt(tau.var.mnar3 / nrow(dat.mnar)) # 95% CI


## IF and MAR-based efficient estimator
IF.mnar.eff1 <- mu.mnar.eff1.1 - mu.mnar.eff1.0 + 
  ifelse(dat.mnar$R == 1 | dat.mnar$S1 == 1, 
         (dat.mnar$A * (dat.mnar$ptwc - mu.mnar.eff1.1) / pi.cc - 
            (1 - dat.mnar$A) * (dat.mnar$ptwc - mu.mnar.eff1.0) / (1 - pi.cc)), 0) / 
  (1 - (1 - gamma.mnar) * (1 - eta.0.mnar.1))

taustar.mnar1 <- mean(IF.mnar.eff1)
taustar.var.mnar1 <- mean((IF.mnar.eff1 - taustar.mnar1)^2)
taustar.mnar1                                                                 # point estimate
taustar.mnar1 + c(-1, 1)*qnorm(0.975) * sqrt(taustar.var.mnar1 / nrow(dat.mnar)) # 95% CI

IF.mnar.eff2 <- mu.mnar.eff2.1 - mu.mnar.eff2.0 + 
  ifelse(dat.mnar$R == 1 | dat.mnar$S2 == 1, 
         (dat.mnar$A * (dat.mnar$ptwc - mu.mnar.eff2.1) / pi.cc - 
            (1 - dat.mnar$A) * (dat.mnar$ptwc - mu.mnar.eff2.0) / (1 - pi.cc)), 0) / 
  (1 - (1 - gamma.mnar) * (1 - eta.0.mnar.2))

taustar.mnar2 <- mean(IF.mnar.eff2)
taustar.var.mnar2 <- mean((IF.mnar.eff2 - taustar.mnar2)^2)
taustar.mnar2                                                                 # point estimate
taustar.mnar2 + c(-1, 1)*qnorm(0.975) * sqrt(taustar.var.mnar2 / nrow(dat.mnar)) # 95% CI

IF.mnar.eff3 <- mu.mnar.eff3.1 - mu.mnar.eff3.0 + 
  ifelse(dat.mnar$R == 1 | dat.mnar$S3 == 1, 
         (dat.mnar$A * (dat.mnar$ptwc - mu.mnar.eff3.1) / pi.cc - 
            (1 - dat.mnar$A) * (dat.mnar$ptwc - mu.mnar.eff3.0) / (1 - pi.cc)), 0) / 
  (1 - (1 - gamma.mnar) * (1 - eta.0.mnar.3))

taustar.mnar3 <- mean(IF.mnar.eff3)
taustar.var.mnar3 <- mean((IF.mnar.eff3 - taustar.mnar3)^2)
taustar.mnar3                                                                 # point estimate
taustar.mnar3 + c(-1, 1)*qnorm(0.975) * sqrt(taustar.var.mnar3 / nrow(dat.mnar)) # 95% CI


## adaptive estimator
diff.var.mnar1 <- mean(((IF.mnar.np1 - IF.mnar.eff1) - (tau.mnar1 - taustar.mnar1))^2)
R.mod.mnar1 <- c(tau.var.mnar1 / nrow(dat.mnar), taustar.var.mnar1 / nrow(dat.mnar) + 
                  max((tau.mnar1 - taustar.mnar1)^2 - diff.var.mnar1 / nrow(dat.mnar), 0))
cross.cov.mnar1 <- mean((IF.mnar.np1 - tau.mnar1) * (IF.mnar.eff1 -  taustar.mnar1))
tauadapt.mnar1 <- ifelse(R.mod.mnar1[1] < R.mod.mnar1[2], tau.mnar1, taustar.mnar1)

diff.var.mnar2 <- mean(((IF.mnar.np2 - IF.mnar.eff2) - (tau.mnar2 - taustar.mnar2))^2)
R.mod.mnar2 <- c(tau.var.mnar2 / nrow(dat.mnar), taustar.var.mnar2 / nrow(dat.mnar) + 
                  max((tau.mnar2 - taustar.mnar2)^2 - diff.var.mnar2 / nrow(dat.mnar), 0))
cross.cov.mnar2 <- mean((IF.mnar.np2 - tau.mnar2) * (IF.mnar.eff2 -  taustar.mnar2))
tauadapt.mnar2 <- ifelse(R.mod.mnar2[1] < R.mod.mnar2[2], tau.mnar2, taustar.mnar2)

diff.var.mnar3 <- mean(((IF.mnar.np3 - IF.mnar.eff3) - (tau.mnar3 - taustar.mnar3))^2)
R.mod.mnar3 <- c(tau.var.mnar3 / nrow(dat.mnar), taustar.var.mnar3 / nrow(dat.mnar) + 
                  max((tau.mnar3 - taustar.mnar3)^2 - diff.var.mnar3 / nrow(dat.mnar), 0))
cross.cov.mnar3 <- mean((IF.mnar.np3 - tau.mnar3) * (IF.mnar.eff3 -  taustar.mnar3))
tauadapt.mnar3 <- ifelse(R.mod.mnar3[1] < R.mod.mnar3[2], tau.mnar3, taustar.mnar3)

## compute CIs for adaptive estimator
set.seed(264)
M <- 5000
beta.grid <- seq(-300, 300, by = 0.01)

Sigma1 <- rbind(c(tau.var.mnar1, cross.cov.mnar1),
                c(cross.cov.mnar1, taustar.var.mnar1))
Z <- rmvnorm(M, sqrt(nrow(dat.mnar)) * c(tau.mnar1, taustar.mnar1), Sigma1)
## compute indicator of the event A_0^{est}
A.g <- 1 * (tau.var.mnar1 < pmax((Z[,2] - Z[,1])^2, 0) + taustar.var.mnar1)
b.beta <- sapply(beta.grid, function(x) {
  mean((Z[,1] - sqrt(nrow(dat.mnar)) * tau.mnar1 <= x) & A.g == 1) +
    mean((Z[,2] - sqrt(nrow(dat.mnar)) * tau.mnar1 <= x) & A.g == 0)
}, simplify = 0)
b.025 <- beta.grid[which.min(abs(b.beta - 0.025))]
b.975 <- beta.grid[which.min(abs(b.beta - 0.975))]
tauadapt.CI.mnar1 <- c(tauadapt.mnar1 - b.975 / sqrt(nrow(dat.mnar)),
                      tauadapt.mnar1 - b.025 / sqrt(nrow(dat.mnar)))

Sigma2 <- rbind(c(tau.var.mnar2, cross.cov.mnar2),
                c(cross.cov.mnar2, taustar.var.mnar2))
Z <- rmvnorm(M, sqrt(nrow(dat.mnar)) * c(tau.mnar2, taustar.mnar2), Sigma2)
## compute indicator of the event A_0^{est}
A.g <- 1 * (tau.var.mnar2 < pmax((Z[,2] - Z[,1])^2, 0) + taustar.var.mnar2)
b.beta <- sapply(beta.grid, function(x) {
  mean((Z[,1] - sqrt(nrow(dat.mnar)) * tau.mnar2 <= x) & A.g == 1) +
    mean((Z[,2] - sqrt(nrow(dat.mnar)) * tau.mnar2 <= x) & A.g == 0)
}, simplify = 0)
b.025 <- beta.grid[which.min(abs(b.beta - 0.025))]
b.975 <- beta.grid[which.min(abs(b.beta - 0.975))]
tauadapt.CI.mnar2 <- c(tauadapt.mnar2 - b.975 / sqrt(nrow(dat.mnar)),
                      tauadapt.mnar2 - b.025 / sqrt(nrow(dat.mnar)))

Sigma3 <- rbind(c(tau.var.mnar3, cross.cov.mnar3),
                c(cross.cov.mnar3, taustar.var.mnar3))
Z <- rmvnorm(M, sqrt(nrow(dat.mnar)) * c(tau.mnar3, taustar.mnar3), Sigma3)
## compute indicator of the event A_0^{est}
A.g <- 1 * (tau.var.mnar3 < pmax((Z[,2] - Z[,1])^2, 0) + taustar.var.mnar3)
b.beta <- sapply(beta.grid, function(x) {
  mean((Z[,1] - sqrt(nrow(dat.mnar)) * tau.mnar3 <= x) & A.g == 1) +
    mean((Z[,2] - sqrt(nrow(dat.mnar)) * tau.mnar3 <= x) & A.g == 0)
}, simplify = 0)
b.025 <- beta.grid[which.min(abs(b.beta - 0.025))]
b.975 <- beta.grid[which.min(abs(b.beta - 0.975))]
tauadapt.CI.mnar3 <- c(tauadapt.mnar3 - b.975 / sqrt(nrow(dat.mnar)),
                      tauadapt.mnar3 - b.025 / sqrt(nrow(dat.mnar)))



# save(dat.cc, dat.mar, dat.mnar,
#      chi.cc, chi.var.cc, chi.mar, chi.var.mar, chi.mnar, chi.var.mnar,
#      tau.mar1, tau.var.mar1, tau.mnar1, tau.var.mnar1,
#      tau.mar2, tau.var.mar2, tau.mnar2, tau.var.mnar2,
#      tau.mar3, tau.var.mar3, tau.mnar3, tau.var.mnar3,
#      taustar.mar1, taustar.var.mar1, taustar.mnar1, taustar.var.mnar1,
#      taustar.mar2, taustar.var.mar2, taustar.mnar2, taustar.var.mnar2,
#      taustar.mar3, taustar.var.mar3, taustar.mnar3, taustar.var.mnar3,
#      diff.var.mar1, R.mod.mar1, cross.cov.mar1, tauadapt.mar1,
#      diff.var.mar2, R.mod.mar2, cross.cov.mar2, tauadapt.mar2,
#      diff.var.mar3, R.mod.mar3, cross.cov.mar3, tauadapt.mar3,
#      diff.var.mnar1, R.mod.mnar1, cross.cov.mnar1, tauadapt.mnar1,
#      diff.var.mnar2, R.mod.mnar2, cross.cov.mnar2, tauadapt.mnar2,
#      diff.var.mnar3, R.mod.mnar3, cross.cov.mnar3, tauadapt.mnar3,
#      tauadapt.CI.mar1, tauadapt.CI.mnar1,
#      tauadapt.CI.mar2, tauadapt.CI.mnar2,
#      tauadapt.CI.mar3, tauadapt.CI.mnar3,
#      file = "res2.Rdata")
# load("res2.Rdata")


### Results summary
par(mar = c(3,8,1.6,1), mfrow = c(3,1))
# par(mar = c(3,4,1.6,1), mfrow = c(1,1))
plot(NULL, ylim = c(-9,-6.0), xlim = c(0.1,0.9), 
     xlab = "", ylab = "Average Treatment Effect", xaxt = 'n')
axis(1, at = c(0.167, 0.333, 0.500, 0.667, 0.833), 
     labels = c("Full Data", "MAR-incomplete", "NP-eff", "MAR-eff", "Adaptive"))
abline(h = chi.cc, col = 'black', lty = 'dashed')
arrows(x0=0.167, y0=chi.cc - qnorm(0.975) * sqrt(chi.var.cc / nrow(dat.cc)), 
       x1=0.167, y1=chi.cc + qnorm(0.975) * sqrt(chi.var.cc / nrow(dat.cc)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'black')
points(0.167, chi.cc, pch = 20)
arrows(x0=0.313, y0=chi.mar - qnorm(0.975) * sqrt(chi.var.mar / nrow(dat.mar)), 
       x1=0.313, y1=chi.mar + qnorm(0.975) * sqrt(chi.var.mar / nrow(dat.mar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'blue')
points(0.313, chi.mar, pch = 20)
arrows(x0=0.353, y0=chi.mnar - qnorm(0.975) * sqrt(chi.var.mnar / nrow(dat.mnar)), 
       x1=0.353, y1=chi.mnar + qnorm(0.975) * sqrt(chi.var.mnar / nrow(dat.mnar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'red')
points(0.353, chi.mnar, pch = 20)
arrows(x0=0.48, y0=tau.mar1 - qnorm(0.975) * sqrt(tau.var.mar1 / nrow(dat.mar)), 
       x1=0.48, y1=tau.mar1 + qnorm(0.975) * sqrt(tau.var.mar1 / nrow(dat.mar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'blue')
points(0.48, tau.mar1, pch = 20)
arrows(x0=0.52, y0=tau.mnar1 - qnorm(0.975) * sqrt(tau.var.mnar1 / nrow(dat.mnar)), 
       x1=0.52, y1=tau.mnar1 + qnorm(0.975) * sqrt(tau.var.mnar1 / nrow(dat.mnar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'red')
points(0.52, tau.mnar1, pch = 20)
arrows(x0=0.647, y0=taustar.mar1 - qnorm(0.975) * sqrt(taustar.var.mar1 / nrow(dat.mar)), 
       x1=0.647, y1=taustar.mar1 + qnorm(0.975) * sqrt(taustar.var.mar1 / nrow(dat.mar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'blue')
points(0.647, taustar.mar1, pch = 20)
arrows(x0=0.687, y0=taustar.mnar1 - qnorm(0.975) * sqrt(taustar.var.mnar1 / nrow(dat.mnar)), 
       x1=0.687, y1=taustar.mnar1 + qnorm(0.975) * sqrt(taustar.var.mnar1 / nrow(dat.mnar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'red')
points(0.687, taustar.mnar1, pch = 20)
arrows(x0=0.813, y0=tauadapt.CI.mar1[1], 
       x1=0.813, y1=tauadapt.CI.mar1[2], 
       code=3, angle=90, length=0.05, lwd=2,col = 'blue')
points(0.813, tauadapt.mar1, pch = 20)
arrows(x0=0.853, y0=tauadapt.CI.mnar1[1], 
       x1=0.853, y1=tauadapt.CI.mnar1[2], 
       code=3, angle=90, length=0.05, lwd=2,col = 'red')
points(0.853, tauadapt.mnar1, pch = 20)
legend("bottomleft", lty = 1, col = c("black", "blue", "red"), legend = c("Full data", "MAR", "MNAR"),
       bty = 'n')
mtext("a) ", padj=13, adj=-0.18, cex=1.0, font=1)

plot(NULL, ylim = c(-9,-6.0), xlim = c(0.1,0.9), 
     xlab = "", ylab = "Average Treatment Effect", xaxt = 'n')
axis(1, at = c(0.167, 0.333, 0.500, 0.667, 0.833), 
     labels = c("Full Data", "MAR-incomplete", "NP-eff", "MAR-eff", "Adaptive"))
abline(h = chi.cc, col = 'black', lty = 'dashed')
arrows(x0=0.167, y0=chi.cc - qnorm(0.975) * sqrt(chi.var.cc / nrow(dat.cc)), 
       x1=0.167, y1=chi.cc + qnorm(0.975) * sqrt(chi.var.cc / nrow(dat.cc)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'black')
points(0.167, chi.cc, pch = 20)
arrows(x0=0.313, y0=chi.mar - qnorm(0.975) * sqrt(chi.var.mar / nrow(dat.mar)), 
       x1=0.313, y1=chi.mar + qnorm(0.975) * sqrt(chi.var.mar / nrow(dat.mar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'blue')
points(0.313, chi.mar, pch = 20)
arrows(x0=0.353, y0=chi.mnar - qnorm(0.975) * sqrt(chi.var.mnar / nrow(dat.mnar)), 
       x1=0.353, y1=chi.mnar + qnorm(0.975) * sqrt(chi.var.mnar / nrow(dat.mnar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'red')
points(0.353, chi.mnar, pch = 20)
arrows(x0=0.48, y0=tau.mar2 - qnorm(0.975) * sqrt(tau.var.mar2 / nrow(dat.mar)), 
       x1=0.48, y1=tau.mar2 + qnorm(0.975) * sqrt(tau.var.mar2 / nrow(dat.mar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'blue')
points(0.48, tau.mar2, pch = 20)
arrows(x0=0.52, y0=tau.mnar2 - qnorm(0.975) * sqrt(tau.var.mnar2 / nrow(dat.mnar)), 
       x1=0.52, y1=tau.mnar2 + qnorm(0.975) * sqrt(tau.var.mnar2 / nrow(dat.mnar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'red')
points(0.52, tau.mnar2, pch = 20)
arrows(x0=0.647, y0=taustar.mar2 - qnorm(0.975) * sqrt(taustar.var.mar2 / nrow(dat.mar)), 
       x1=0.647, y1=taustar.mar2 + qnorm(0.975) * sqrt(taustar.var.mar2 / nrow(dat.mar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'blue')
points(0.647, taustar.mar2, pch = 20)
arrows(x0=0.687, y0=taustar.mnar2 - qnorm(0.975) * sqrt(taustar.var.mnar2 / nrow(dat.mnar)), 
       x1=0.687, y1=taustar.mnar2 + qnorm(0.975) * sqrt(taustar.var.mnar2 / nrow(dat.mnar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'red')
points(0.687, taustar.mnar2, pch = 20)
arrows(x0=0.813, y0=tauadapt.CI.mar2[1], 
       x1=0.813, y1=tauadapt.CI.mar2[2], 
       code=3, angle=90, length=0.05, lwd=2,col = 'blue')
points(0.813, tauadapt.mar2, pch = 20)
arrows(x0=0.853, y0=tauadapt.CI.mnar2[1], 
       x1=0.853, y1=tauadapt.CI.mnar2[2], 
       code=3, angle=90, length=0.05, lwd=2,col = 'red')
points(0.853, tauadapt.mnar2, pch = 20)
legend("bottomleft", lty = 1, col = c("black", "blue", "red"), legend = c("Full data", "MAR", "MNAR"),
       bty = 'n')
mtext("b) ", padj=13, adj=-0.18, cex=1.0, font=1)

plot(NULL, ylim = c(-9,-6.0), xlim = c(0.1,0.9), 
     xlab = "", ylab = "Average Treatment Effect", xaxt = 'n')
axis(1, at = c(0.167, 0.333, 0.500, 0.667, 0.833), 
     labels = c("Full Data", "MAR-incomplete", "NP-eff", "MAR-eff", "Adaptive"))
abline(h = chi.cc, col = 'black', lty = 'dashed')
arrows(x0=0.167, y0=chi.cc - qnorm(0.975) * sqrt(chi.var.cc / nrow(dat.cc)), 
       x1=0.167, y1=chi.cc + qnorm(0.975) * sqrt(chi.var.cc / nrow(dat.cc)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'black')
points(0.167, chi.cc, pch = 20)
arrows(x0=0.313, y0=chi.mar - qnorm(0.975) * sqrt(chi.var.mar / nrow(dat.mar)), 
       x1=0.313, y1=chi.mar + qnorm(0.975) * sqrt(chi.var.mar / nrow(dat.mar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'blue')
points(0.313, chi.mar, pch = 20)
arrows(x0=0.353, y0=chi.mnar - qnorm(0.975) * sqrt(chi.var.mnar / nrow(dat.mnar)), 
       x1=0.353, y1=chi.mnar + qnorm(0.975) * sqrt(chi.var.mnar / nrow(dat.mnar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'red')
points(0.353, chi.mnar, pch = 20)
arrows(x0=0.48, y0=tau.mar3 - qnorm(0.975) * sqrt(tau.var.mar3 / nrow(dat.mar)), 
       x1=0.48, y1=tau.mar3 + qnorm(0.975) * sqrt(tau.var.mar3 / nrow(dat.mar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'blue')
points(0.48, tau.mar3, pch = 20)
arrows(x0=0.52, y0=tau.mnar3 - qnorm(0.975) * sqrt(tau.var.mnar3 / nrow(dat.mnar)), 
       x1=0.52, y1=tau.mnar3 + qnorm(0.975) * sqrt(tau.var.mnar3 / nrow(dat.mnar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'red')
points(0.52, tau.mnar3, pch = 20)
arrows(x0=0.647, y0=taustar.mar3 - qnorm(0.975) * sqrt(taustar.var.mar3 / nrow(dat.mar)), 
       x1=0.647, y1=taustar.mar3 + qnorm(0.975) * sqrt(taustar.var.mar3 / nrow(dat.mar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'blue')
points(0.647, taustar.mar3, pch = 20)
arrows(x0=0.687, y0=taustar.mnar3 - qnorm(0.975) * sqrt(taustar.var.mnar3 / nrow(dat.mnar)), 
       x1=0.687, y1=taustar.mnar3 + qnorm(0.975) * sqrt(taustar.var.mnar3 / nrow(dat.mnar)), 
       code=3, angle=90, length=0.05, lwd=2,col = 'red')
points(0.687, taustar.mnar3, pch = 20)
arrows(x0=0.813, y0=tauadapt.CI.mar3[1], 
       x1=0.813, y1=tauadapt.CI.mar3[2], 
       code=3, angle=90, length=0.05, lwd=2,col = 'blue')
points(0.813, tauadapt.mar3, pch = 20)
arrows(x0=0.853, y0=tauadapt.CI.mnar3[1], 
       x1=0.853, y1=tauadapt.CI.mnar3[2], 
       code=3, angle=90, length=0.05, lwd=2,col = 'red')
points(0.853, tauadapt.mnar3, pch = 20)
legend("bottomleft", lty = 1, col = c("black", "blue", "red"), legend = c("Full data", "MAR", "MNAR"),
       bty = 'n')
mtext("c) ", padj=13, adj=-0.18, cex=1.0, font=1)

## save as portrait, 6 x 11 inches
