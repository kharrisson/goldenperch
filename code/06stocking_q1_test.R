library(dplyr)
library(rstanarm)
data <- read.csv("./data/stocking_q1.csv")
data <- data |>
  mutate(nstocked_std = (n_stocked - mean(n_stocked)) / sd(n_stocked))

# The number of stocked fish was z-scaled prior to analysis.

mod_lm <- stan_lm(q1 ~ nstocked_std, data, prior = R2(location = 0.5, what = "median"))

# There was evidence of a negative association between Q1 and the number of stocked fish
#.  (mean scaled estimate equal to X [90% CrI x-y])
# The mean estimate was
estimate <- posterior_interval(mod_lm, prob = 0.01)
round(mean(c(estimate["nstocked_std", 1], estimate["nstocked_std", 2])), 2)
# (90% credible interval: )
cri <- posterior_interval(mod_lm, prob = 0.9)
round(c(cri["nstocked_std", 1], cri["nstocked_std", 2]), 2)
