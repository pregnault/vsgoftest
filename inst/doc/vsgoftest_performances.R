#### header ####
## vsgoftest_performances.R
## Code chunks for generating outputs presented in Section "Performances of Vasicek-Song tests" of 
## J Lequesne, P Regnault (2018). Package vsgoftest for R: goodness-of-fit tests based on Kullback-Leibler divergence.

#### List of required packages ####

library(vsgoftest)
library(dbEmpLikeGOF)
library(ggplot2)
library(dplyr)
library(knitr)
library(microbenchmark)
library(goftest)

#### vsgoftest versus dbEmpLikeGOF ####

## For normal distribution
set.seed(1)
bm <- microbenchmark(vs.test = vs.test(x = (sample <- rnorm(n = 50)), densfun = 'dnorm', 
                                       simulate.p.value = TRUE, B = 1000, delta = -1/6),
                     dbEmpLikeGOF = dbEmpLikeGOF(x = sample, testcall = "normal", 
                                                 pvl.Table = FALSE, num.mc = 1000, vrb = FALSE),
                     times = 100L)
bm %>%
  rename(`function` = expr) %>%
  group_by(`function`) %>%
  summarize(min = min(time/10^6),
            `1stQ` = quantile(time/10^6,0.25),
            median = median(time/10^6),
            mean = mean(time/10^6),
            `3rdQ` = quantile(time/10^6,0.75),
            max = max(time/10^6),
            sd = sd(time/10^6)) %>%
  kable()

# Violin plots for comparing computation times of vs.test and dbEmpLikeGOF
# postscript(file = '../../papier/PKGFigCompTimeNorm.ps')
bm %>% 
  ggplot(mapping = aes(x = expr, y = time/10^6)) +
  geom_violin() +
  coord_flip() +
  stat_summary(fun.y = "median", geom = "errorbar", 
               mapping = aes(ymax = ..y.., ymin = ..y..),
               linetype = "dashed", col = 'red', show.legend = TRUE) +
  labs(x = '', y = 'Computation time (ms)') +
  theme(text = element_text(size = 18))
# dev.off()

## For uniform distribution
set.seed(1)
bm <- microbenchmark(vs.test = vs.test(x = (sample <- runif(n = 50, min = 1, max = 3)), densfun = 'dunif', 
                                       simulate.p.value = TRUE, B = 1000, delta = -1/6),
                     dbEmpLikeGOF = dbEmpLikeGOF(x = sample, testcall = "uniform",  
                                                 pvl.Table = FALSE, num.mc = 1000, vrb = FALSE),
                     times = 100L)
bm %>%
  rename(`function` = expr) %>%
  group_by(`function`) %>%
  summarize(min = min(time/10^6),
            `1stQ` = quantile(time/10^6,0.25),
            median = median(time/10^6),
            mean = mean(time/10^6),
            `3rdQ` = quantile(time/10^6,0.75),
            max = max(time/10^6),
            sd = sd(time/10^6)) %>%
  kable()

# Violin plots for comparing computation times of vs.test and dbEmpLikeGOF
# postscript(file = '../../papier/PKGFigCompTimeUnif.ps')
bm %>% 
  ggplot(mapping = aes(x = expr, y = time/10^6)) +
  geom_violin() +
  coord_flip() +
  stat_summary(fun.y = "median", geom = "errorbar", 
               mapping = aes(ymax = ..y.., ymin = ..y..),
               linetype = "dashed", col = 'red', show.legend = TRUE) +
  labs(x = '', y = 'Computation time (ms)') +
  theme(text = element_text(size = 18))
# dev.off()


#### Power comparisons ####

lstn <- c(20, 30, 50, 100) #List of samples sizes
N <- 10000 #Number of replicates for MC simulations


## Pareto against log-normal

# auxiliary function for simulating a sample, applying GOF tests and gathering their pvalues
pop <- function(n, m = 0, s) {
  ech <- 1 + rlnorm(n, meanlog = m, sdlog = s)
  pvs <- vs.test(x = ech, 
                 densfun = 'dpareto', param = c(1/s, 1),
                 B = 1000)$p.value
  pks <- ks.test(x = ech, y = 'ppareto', mu = 1/s, c = 1)$p.value
  pad <- ad.test(x = ech, null = 'ppareto', mu = 1/s, c = 1)$p.value
  pcvm <- cvm.test(x = ech, null = 'ppareto', mu = 1/s, c = 1)$p.value
  return(c(pvs, pks, pad, pcvm))
}

#Pareto with c = 1, mu = 1 against shifted log-normal with meanlog = 0, sdlog = 1
mu <- 1
powers1 <- matrix(0, nrow = 4, ncol = length(lstn))
set.seed(54)
for (i in seq_along(lstn)) {
  res.pow <- replicate(n = N, expr = pop(lstn[i], 0, 1/mu))
  powers1[i,] <- apply(res.pow < 0.05, 1, mean)
}

#Pareto with c = 1, mu = 0.8 against shifted log-normal with meanlog = 0, sdlog = 1.25
mu <- 4/5
powers2 <- matrix(0, nrow = 4, ncol = length(lstn))
set.seed(32)
for (i in seq_along(lstn)) {
  res.pow <- replicate(n = N, expr = pop(lstn[i], 0, 1/mu))
  powers2[i,] <- apply(res.pow < 0.05, 1, mean)
}


## Exponential vs Weibull

#Auxialiary function
pop <- function(n, sh) {
  ech <- rweibull(n, shape = sh, scale = 1)
  pvs <- vs.test(x = ech, 
                 densfun = 'dexp', param = 1,
                 B = 1000)$p.value
  pks <- ks.test(x = ech, y = 'pexp', rate = 1)$p.value
  pad <- ad.test(x = ech, null = 'pexp', rate = 1)$p.value
  pcvm <- cvm.test(x = ech, null = 'pexp', rate = 1)$p.value
  return(c(pvs, pks, pad, pcvm))
}

# Exp 1/2 vs Weibull(1.2, 2)
powers3 <- matrix(0, nrow = 4, ncol = length(lstn))
set.seed(12)
for (i in seq_along(lstn)) {
  res.pow <- replicate(n = N, expr = pop(lstn[i], 1.2))
  powers3[i,] <- apply(res.pow < 0.05, 1, mean)
}

# Exp 1/2 vs Weibull(1.3, 2)
powers4 <- matrix(0, nrow = 4, ncol = length(lstn))
set.seed(23)
for (i in seq_along(lstn)) {
  res.pow <- replicate(n = N, expr = pop(lstn[i], 1.3))
  powers4[i,] <- apply(res.pow < 0.05, 1, mean)
}
