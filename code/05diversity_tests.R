
# set working directory
setwd("~/Dropbox/gp_2022")

library(dplyr)

div <- read.csv("./manuscript/supp_files/individual_data_plus_cov.csv")

boxplot(PHt.y~inferred_stocked_conservative.1, div)
anova1 <- aov(PHt.y~inferred_stocked_conservative.1, div)
summary(anova1)
tukey <- TukeyHSD(anova1)

ggplot(div, aes(x = inferred_stocked_conservative.1, y = PHt.y, fill = inferred_stocked_conservative.1)) +
  geom_boxplot() +
  geom_jitter(shape = 15,
              color = "steelblue",
              position = position_jitter(0.21)) +
  theme_classic()

ggplot(div, aes(x = pop.x, y = PHt.y, fill = pop.code)) +
  geom_boxplot() 
  #+
  #geom_jitter(shape = 15,
             # color = "steelblue",
             # position = position_jitter(0.21)) +
  theme_classic()


mean_PHt = mean(div$PHt.y, na.rm = TRUE),
sd_PHt = sd(div$PHt.y, na.rm = TRUE)

div %>%
  group_by(inferred_stocked_conservative.1) %>%
  summarise_at(vars(PHt.y), list(name = mean))



div_known_origin <- div[ which(div$inferred_stocked_conservative.1=='YES'
                         | div$inferred_stocked_conservative.1 == "NO"), ]

PHt_stocked_plot <- div_known_origin %>%
  ggplot(aes(x = inferred_stocked_conservative.1, y = PHt.y)) +
  geom_boxplot(fill="gray", outlier.shape=1) +
  labs(title="",x="stocking origin", y = "PHt") +
  geom_jitter(
    shape = 16,
    color = "black",
    position = position_jitter(0.1)
  ) +
  theme_classic()
ggsave(
  PHt_stocked_plot,
  filename = "figures/figure_PHt_stocked.png",
  device = png,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

div_known_origin %>%
  group_by(inferred_stocked_conservative.1) %>%
  summarise_at(vars(PHt.y), list(name = mean))

ttest <- t.test(PHt.y ~ inferred_stocked_conservative.1, data = div_known_origin)

lme_test <- lmerTest::lmer(
  PHt.y ~ inferred_stocked_conservative.1 + 
    (1 | pop.x), 
  data = div_known_origin
)

