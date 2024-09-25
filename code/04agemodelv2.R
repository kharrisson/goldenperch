#set working directory
setwd("~/Dropbox/gp_2022")

#load packages
library(qs)
library(betareg)
library(rstanarm)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)
library(patchwork)

gp_age_data <- read.csv("./data/structure_output_aged_fish.csv", header =TRUE)
gp_cov <- read.csv("./data/individual_data_plus_cov.csv", header =TRUE)
  
#make a new dataframe
gp_age_data_clean <- data.frame (id  = gp_age_data$id,
                             age = gp_age_data$oto_age,
                             q1 = gp_age_data$q1,
                             sibling_group = as.factor(gp_age_data$sibling_group),
                             pop = as.factor(gp_age_data$pop))

#specify number of chains
num_chains <- 4
num_iter <- 10000

#fit poisson regression models
gp_data_clean_filled <- gp_age_data_clean %>%
  mutate(age_std = scale(age)  ## standardise so it has a mean of zero and SD = 1.
  )

#each individual needs to be its only thing if it's not part of a sibling group - i.e. not all
#unrelated, will add lots of parameters = bad but oh well. 
idx <- is.na(gp_data_clean_filled$sibling_group)
gp_data_clean_filled$sibling_group <- as.numeric(gp_data_clean_filled$sibling_group)
max_val <- max(gp_data_clean_filled$sibling_group[!idx], na.rm = TRUE)
gp_data_clean_filled$sibling_group[idx] <- max_val + seq_len(sum(idx))
gp_data_clean_filled$sibling_group <- factor(gp_data_clean_filled$sibling_group)

## CHANGE THIS SO YOU GET WARNINGS AS WARNING RATHER THAN CHANGING
##   THEM TO ERRORS
options(warn = 0)

# add a small amount to the q value
eps <- min(
  gp_data_clean_filled %>% filter(q1 > 0) %>% pull(q1),
  na.rm = TRUE
)
gp_data_clean_filled <- gp_data_clean_filled %>%
  mutate(
    q1 = ifelse(q1 == 0, q1 + (eps / 2.), q1),
    q1 = ifelse(q1 == 1, q1 - (eps / 2.), q1)
  )

mod_age <- stan_glmer(
  q1 ~ age_std +  
    (1 | pop) +
    (1 | sibling_group),
  data = gp_data_clean_filled,
  family = mgcv::betar,
  prior = NULL,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 20
  ),
  chains = 4,#num_chains,
  iter = 10000,#num_iter,
  cores = 4
)


summary(mod_age)
plot(mod_age, prob=c(0.8), prob_outer = 0.95, pars = c("(Intercept)", "age_std"))                
posterior_interval(mod_age, prob=c(0.8), pars = c("(Intercept)", "age_std"))


p1 <- plot(mod_stocked_or_not, prob=c(0.8), prob_outer = 0.95,
           pars = c("(Intercept)", "inferred_stocked.yYES"))  



##effects package (to plot the relationship)
##change to chronological (years)
ages_to_predict <- seq(-5, 26, by = 1)
predicted_q1 <- posterior_predict(
  mod_age, 
  newdata = data.frame(age_std = ages_to_predict),
  re.form = NA
)
effects <- data.frame(
  value = apply(predicted_q1, 2, median),
  lower1 = apply(predicted_q1, 2, quantile, probs = 0.25),
  upper1 = apply(predicted_q1, 2, quantile, probs = 0.75),
  lower2 = apply(predicted_q1, 2, quantile, probs = 0.1),
  upper2 = apply(predicted_q1, 2, quantile, probs = 0.9)
)

raw_data_to_plot <- gp_data_clean_filled %>%
  mutate(birth_year = 2018L - age - 1L)

effects <- effects %>%
  mutate(
    age = ages_to_predict,
    birth_year = 2018L - age - 1L,
    speculative = ifelse(birth_year >= 2016, value, NA),
    value = ifelse(birth_year <= 2016, value, NA)
  )
# speculative_plot <- effects %>%
#   ggplot(
#     aes(x = birth_year, y = value, ymin = lower1, ymax = upper1)
#   ) +
#   geom_line() +
#   geom_line(aes(y = speculative), linetype = "dashed") +
#   geom_ribbon(alpha = 0.3) +
#   geom_ribbon(aes(ymin = lower2, ymax = upper2), alpha = 0.15) +
#   annotate(
#     "point",
#     x = raw_data_to_plot$birth_year, 
#     y = raw_data_to_plot$q1,
#     col = "grey25"
#   ) +
#   xlab("Birth year") +
#   ylab("Admixture coefficient (q1)") +
#   scale_x_continuous(
#     breaks = seq(min(effects$birth_year), max(effects$birth_year), by = 5)
#   )

boring_katherine_plot <- effects %>%
  filter(birth_year <= 2016) %>%
  ggplot(
    aes(x = birth_year, y = value, ymin = lower1, ymax = upper1)
  ) +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  geom_ribbon(aes(ymin = lower2, ymax = upper2), alpha = 0.15) +
  annotate(
    "point",
    x = raw_data_to_plot$birth_year, 
    y = raw_data_to_plot$q1,
    col = "grey25"
  ) +
  xlab("Birth year") +
  ylab("Admixture coefficient (q1)") +
  scale_x_continuous(
    breaks = seq(min(effects$birth_year), 2016, by = 5)
  )


gp_data_stocked <- gp_data_clean_filled %>%
  left_join(gp_cov %>% select(id, inferred_stocked.y)) %>%
  filter(inferred_stocked.y %in% c("NO", "YES"))

gp_data_stocked_yes <- gp_data_stocked %>%
  filter(inferred_stocked.y == "YES")

gp_data_stocked_no <- gp_data_stocked %>%
  filter(inferred_stocked.y == "NO")


mod_stocked_or_not <- stan_glmer(
  q1 ~ inferred_stocked.y +  
    (1 | pop) +
    (1 | sibling_group),
  data = gp_data_stocked,
  family = mgcv::betar,
  prior = NULL,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 20
  ),
  chains = 4,#num_chains,
  iter = 10000,#num_iter,
  cores = 4
)


summary(mod_stocked_or_not)

p1 <- plot(mod_stocked_or_not, prob=c(0.8), prob_outer = 0.95,
                                pars = c("(Intercept)", "inferred_stocked.yYES"))                
p1 + ggplot2::ggtitle("Posterior medians \n with 80% and 95% intervals")
posterior_interval(mod_stocked_or_not, prob=0.8, pars = c("(Intercept)", "inferred_stocked.yYES"))

stocked_predictions <- posterior_predict(
  mod_stocked_or_not, 
  newdata = data.frame(inferred_stocked.y = c("YES", "NO")),
  re.form = NA
)
effects <- data.frame(
  value = apply(stocked_predictions, 2, median),
  lower1 = apply(stocked_predictions, 2, quantile, probs = 0.25),
  upper1 = apply(stocked_predictions, 2, quantile, probs = 0.75),
  lower2 = apply(stocked_predictions, 2, quantile, probs = 0.1),
  upper2 = apply(stocked_predictions, 2, quantile, probs = 0.9)
)
# gp_data_stocked

effects <- effects %>%
  mutate(stocked = c("YES", "NO"))
effects %>% 
  ggplot(aes(x = stocked, y = value, ymin = lower1, ymax = upper1)) +
  geom_point(size = 2) + 
  geom_errorbar(width = 0.25) +
  geom_errorbar(aes(ymin = lower2, ymax = upper2), width = 0.1) +
  annotate(
    "point",
    x = gp_data_stocked$inferred_stocked.y, 
    y = gp_data_stocked$q1,
    col = gp_data_stocked$sibling_group
  ) +
  xlab("") +
  ylab("")

# Build dataset with different distributions

library(ggplot2)
library(dplyr)
library(hrbrthemes)

# Represent it
p <- gp_data_stocked  %>%
  ggplot( aes(x=q1, fill=inferred_stocked.y)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="")

# look at raw data to work out proportion of stocked/not stocked
#    when q1 >= 0.375
stocked_q1_gt375 <- gp_data_stocked |>
  filter(q1 > 0.375) |>
  group_by(inferred_stocked.y) |>
  summarise(count = n()) |>
  mutate(
    total = sum(count),
    proportion = count / total
  )

# now use the model to do something
proportion_modelled_q1_lt75 <- stocked_predictions |>
  apply(2, \(x) sum(x <= 0.75) / length(x))
proportion_modelled_q1_lt75 <- proportion_modelled_q1_lt75 / sum(proportion_modelled_q1_lt75)
names(proportion_modelled_q1_lt75) <- c("YES", "NO")

proportion_modelled_q1_gt75 <- stocked_predictions |>
  apply(2, \(x) sum(x > 0.75) / length(x))
proportion_modelled_q1_gt75 <- proportion_modelled_q1_gt75 / sum(proportion_modelled_q1_gt75)
names(proportion_modelled_q1_gt75) <- c("YES", "NO")
