library(haven)
library(dplyr)

select <- dplyr::select

#---- DATASET STRUCTURE-----
# exposure (X): equal to 1 if a uniformly distributed random number between 0 and 1 is <0.4
# confounder (C): rnormal()*0.75
# mediator (M): 2*X + 1.6*C + rnormal()*0.75
# outcome (Y): Y = 1.5*X + 1*M + .8*C + rnormal()*2

# reading datasets generated from stata
d <- read_dta("data.dta")

# make id column
d$id <- c(1:10000)

# number of bootstraps
nBoot <- 1000

#----CONTROLLED DIRECT EFFECT-----
# vector to save all estimates from bootstrap
cde_est <- c()

for(i in 1:nBoot){
  # select ids
  ids <- sample(d$id, replace = TRUE)
  
  # make bootstrapped sample
  d_boot <- d[ids, ]
  
  # make dataset copies so that each row is repeated 3 times
  d_boot <- d_boot %>%
    slice(rep(1:n(), each=3)) %>%
    mutate(copy = rep(1:3, 10000))
  
  # copy 2: set x = 0, m = 0, y = missing
  d_boot <- d_boot %>%
    mutate(x = ifelse(copy == 2, 0, x),
           m = ifelse(copy == 2, 0, m),
           y = ifelse(copy == 2, NA, y))
  
  # copy 3: set x = 1, m = 0, y = missing
  d_boot <- d_boot %>%
    mutate(x = ifelse(copy == 3, 1, x),
           m = ifelse(copy == 3, 0, m),
           y = ifelse(copy == 3, NA, y))
  
  # estimate model predicting the outcome as a function of exposure, mediator, 
  # the interaction of the exposure and mediator, and the mediator-outcome 
  # confounder (C), using only the real observations
  mod_cde <- lm(y ~ x + m + m:x + c, data = d_boot, subset = (copy == 1))
  
  # counterfactual value of Y when X = 0 and M = 0
  cf_y_x0_m0 <- predict(mod_cde, newdata = d_boot[d_boot$copy == 2, ])
  
  # counterfactual value of Y when X = 1 and M = 1
  cf_y_x1_m0 <- predict(mod_cde, newdata = d_boot[d_boot$copy == 3, ])
  
  # CDE of X on Y
  cde_est[i] <- mean(cf_y_x1_m0) - mean(cf_y_x0_m0)
}

# mean and 95% CI of bootstrapped CDE estimates
mean(cde_est)
quantile(cde_est, probs=c(0.025,0.975))

#----NATURAL INDIRECT EFFECT-----
# vector to save all estimates from bootstrap
nie_est <- c()

# vector to save % variance explained
pve <- c()

for(i in 1:nBoot){
  # select ids
  ids <- sample(d$id, replace = TRUE)
  
  # make bootstrapped sample
  d_boot <- d[ids, ]
  
  # make dataset copies
  d_boot <- d_boot %>%
    slice(rep(1:n(), each=3)) %>%
    mutate(copy = rep(1:3, 10000))
  
  # copy 2: set x = 0, m = missing, y = missing
  d_boot <- d_boot %>%
    mutate(x = ifelse(copy == 2, 0, x),
           m = ifelse(copy == 2, NA, m),
           y = ifelse(copy == 2, NA, y))
  
  # copy 3: set x = 1, m = missing, y = missing
  d_boot <- d_boot %>%
    mutate(x = ifelse(copy == 3, 1, x),
           m = ifelse(copy == 3, NA, m),
           y = ifelse(copy == 3, NA, y))
  
  # estimate model predicting the mediator as a function of X and C, 
  # using only the real observations
  mod_m <- lm(m ~ x + c, data = d_boot, subset = (copy == 1))
  
  # counterfactual value of M when X = 0
  cf_m_x0 <- predict(mod_m, newdata = d_boot[d_boot$copy == 2, ])
  
  # counterfactual value of M when X = 1
  cf_m_x1 <- predict(mod_m, newdata = d_boot[d_boot$copy == 3, ])
  
  # estimate model predicting the outcome as a function of exposure,
  # mediator, the interaction of the exposure and mediator, and the 
  # mediator-outcome confounder, using only the real observations
  mod_y <- lm(y ~ x + m + x:m + c, data = d_boot, subset = (copy == 1))
  
  # set M to cf_m_x0 and set X to 1 in copy = 2
  d_boot$m[d_boot$copy == 2] <- cf_m_x0 
  d_boot$x[d_boot$copy == 2] <- 1
  
  # set M to cf_m_x1 in copy = 3 (note X is already set to 1)
  d_boot$m[d_boot$copy == 3] <- cf_m_x1
  
  # counterfactual value of Y when X=1 and M=the value it would take if X=0. 
  cf_y_x1_cf_m_x0 <- predict(mod_y, newdata = d_boot[d_boot$copy == 2,])
  
  # counterfactual value of Y when X=1 and M=the value it would take if X=1.
  cf_y_x1_cf_m_x1 <- predict(mod_y, newdata = d_boot[d_boot$copy == 3,])
  
  # natural indirect effect of X on Y mediated by M
  nie_est[i] <- mean(cf_y_x1_cf_m_x1) - mean(cf_y_x1_cf_m_x0)
  
  # estimate total effect to calculate % variance explained
  mod_tot <- lm(y ~ x, data = d_boot, subset = (copy == 1))
  tot_eff <- mod_tot$coef['x']
  
  # % variance explained = natural indirect effect/total effect
  pve[i] <- nie_est[i]/tot_eff
  
}

# mean and 95% CI of bootstrapped NIE estimates
mean(nie_est)
quantile(nie_est, probs=c(0.025,0.975))

# mean and 95% CI percent variance explained
mean(pve)
quantile(pve, probs=c(0.025,0.975))

#----NATURAL DIRECT EFFECT-----
# vector to save all estimates
nde_est <- c()

for(i in 1:nBoot){
  # select ids
  ids <- sample(d$id, replace = TRUE)
  
  # make bootstrapped sample
  d_boot <- d[ids, ]
  
  # make dataset copies
  d_boot <- d_boot %>%
    slice(rep(1:n(), each=3)) %>%
    mutate(copy = rep(1:3, 10000))
  
  # copy 2: set x = 0, m = missing, y = missing
  d_boot <- d_boot %>%
    mutate(x = ifelse(copy == 2, 0, x),
           m = ifelse(copy == 2, NA, m),
           y = ifelse(copy == 2, NA, y))
  
  # copy 3: set x = 1, m = missing, y = missing
  d_boot <- d_boot %>%
    mutate(x = ifelse(copy == 3, 1, x),
           m = ifelse(copy == 3, NA, m),
           y = ifelse(copy == 3, NA, y))
  
  # estimate model predicting the mediator as a function of X and C, 
  # using only the real observations
  mod_m <- lm(m ~ x + c, data = d_boot, subset = (copy == 1))
  
  # counterfactual value of M when X = 0
  cf_m_x0 <- predict(mod_m, newdata = d_boot[d_boot$copy == 2,])
  
  # estimate model predicting the outcome as a function of exposure,
  # mediator, the interaction of the exposure and mediator, and the 
  # mediator-outcome confounder, using only the real observations
  mod_y <- lm(y ~ x + m + x:m + c, data = d_boot, subset = (copy == 1))
  
  # set value of M in both copy 2 and 3 to cf_m_x0
  d_boot$m[d_boot$copy == 2] <- cf_m_x0  
  d_boot$m[d_boot$copy == 3] <- cf_m_x0 
  
  # counterfactual value of Y setting X to 0 and M to the value it would take if 
  # X were set to 0.
  cf_y_x0_cf_m_x0 <- predict(mod_y, newdata = d_boot[d_boot$copy == 2,])
  
  # counterfactual value of Y setting X to 1 and M to the value it would take if 
  # X were set to 0.
  cf_y_x1_cf_m_x0 <- predict(mod_y, newdata = d_boot[d_boot$copy == 3,])
  
  # natural direct effect
  nde_est[i] <- mean(cf_y_x1_cf_m_x0) - mean(cf_y_x0_cf_m_x0)
}

# mean and 95% CI of bootstrapped NDE estimates
mean(nde_est)
quantile(nde_est, probs=c(0.025,0.975))

#---- USING R PACKAGE (mediation)------
# ACME = NIE (i used x = 1 above, so ACME (treated) should match up w/ NIE)
# ADE = NDE (i used x = 0 above, so ADE (control) should match up w/ NDE)

library(mediation)

model.m <- lm(m ~ x + c, data = d)
model.y <- lm(y ~ x + m + x:m + c, data = d)

out <- mediate(model.m, model.y, treat = 'x', mediator = 'm', boot = T)
summary(out)


