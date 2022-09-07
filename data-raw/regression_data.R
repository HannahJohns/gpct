# Continuous data with varying strengths

rm(list=ls())
library(tidyverse)
set.seed(-9225552)

n_per_group <- 1000

data.frame(
  coef = rep(1:5,each=n_per_group)
) %>%
  group_by(coef) %>%
    mutate(i=1:n(),
           x=seq(0,1,length.out=n()),
           y=coef*x+rnorm(n(),mean=0,sd=0.1)) %>%
  ungroup() -> interaction_regression

data.frame(
  coef = rep(1:5,each=n_per_group)
) %>%
  group_by(coef) %>%
  mutate(i=1:n(),
         x=seq(0,1,length.out=n()),
         y=3*x+coef+rnorm(n(),mean=0,sd=0.1)) %>%
  ungroup() -> multiple_regression


usethis::use_data(interaction_regression, overwrite = TRUE)
usethis::use_data(multiple_regression, overwrite = TRUE)


# ggplot(multiple_regression,aes(x=x,y=y,color=factor(coef))) + geom_point()
# ggplot(interaction_regression,aes(x=x,y=y,color=factor(coef))) + geom_point()
