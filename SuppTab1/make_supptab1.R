# generate supplementary table 1

library(xtable)
library(dplyr)
library(magrittr)

pi1Vec <- seq(from=0.02, to=0.08, by=0.02)
muVec <- c(3)
rhoVec <- c(0.1, 0.5, 0.8)
paramDF <- expand.grid(pi1Vec, muVec, rhoVec) %>%
  set_colnames(c("pi1", "mu", "rho")) %>%
  mutate(pi3 = pi1^2 * 3) %>%
  mutate(pi0 = 1 - 2*pi1 - pi3) %>%
  select(pi1, pi0, pi3, rho, mu) %>%
  mutate(Cov  = NA, Var=NA)

for (row_it in 1:nrow(paramDF)) {
 
  rho <- paramDF$rho[row_it]
  pi0 <- paramDF$pi0[row_it] 
  pi1 <- paramDF$pi1[row_it]
  pi2 <- pi1
  pi3 <- paramDF$pi3[row_it]
  mu31 <- paramDF$mu[row_it]
  mu32 <- mu31
  mu1 <- mu31
  mu2 <- mu31
  expCov <- rho + pi3 * mu31 * mu32 - (pi1*mu1 + pi3*mu31) * (pi2*mu2 + pi3*mu32)
  expVar <- pi1 * (mu1^2 + 1) + pi3 * (mu31^2 + 1) + pi0 * 1 + pi2 * 1 - (pi1 * mu1 + pi3 * mu31)^2
  paramDF$Cov[row_it] <- expCov
  paramDF$Var[row_it] <- expVar
}
paramDF$Corr <- paramDF$Cov / paramDF$Var

#print(xtable(paramDF, digits=c(0, 2, 2, 4, 2, 0, 2)), include.rownames=F)
print(xtable(paramDF %>% select(-pi0, -mu), digits=c(0, 2,  4, 2, 2, 2, 2)), include.rownames=F)








