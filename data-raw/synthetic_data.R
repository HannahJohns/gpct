# Continuous data with varying strengths

rm(list=ls())
library(tidyverse)
# library(gpct)

n_per_group <- 500
epsilon <- 0.5
coefs <- c(0,1,2)


{
set.seed(-78975072)


data.frame(
  coef = rep(coefs,each=n_per_group),
  type="varying intercept"
) %>%
  group_by(coef) %>%
  mutate(i=1:n(),
         x=seq(0,1,length.out=n()),
         y=x+coef+rnorm(n(),mean=0,sd=epsilon)) %>%
  ungroup() %>%
  mutate(
    x=x,
    formula = sprintf("Y~X+%d+N(0,%0.1f)",coef,epsilon)
    )-> varying_intercept

data.frame(
  coef = rep(coefs,each=n_per_group)
) %>%
  group_by(coef) %>%
  mutate(type="varying slope",
         i=1:n(),
         x=seq(0,1,length.out=n()),
         y=coef*x,
         formula = sprintf("Y~%dX+N(0,%0.1f)",coef,epsilon)
         ) %>%
  ungroup() %>% rbind() %>%
  rbind(.,
    data.frame(coef=rep(10,n_per_group),
               type="varying slope",
               i = 1:n_per_group,
               x = seq(0,1,length.out=n_per_group)
               ) %>% mutate(y= 10*x*(1-x),
                            formula = sprintf("Y~10X*(1-X)+N(0,%0.1f)",epsilon)
                            )
  ) %>%
  mutate(y=y+rnorm(n(),mean=0,sd=epsilon)) -> varying_slope



data.frame(
  epsilon = rep(c(0.5,1,2),each=n_per_group),
  type="varying spread"
) %>%
  group_by(epsilon) %>%
  mutate(i=1:n(),
         x=seq(0,1,length.out=n()),
         y=x+rnorm(n(),mean=0,sd=epsilon)) %>%
  ungroup() %>%
  mutate(
    x=x,
    formula = sprintf("Y~X+N(0,%0.1f)",epsilon)
  )-> varying_spread

}

formula_levels <- unique(c(unique(varying_intercept$formula),unique(varying_slope$formula),unique(varying_spread$formula)))

formula_labels <- gsub("\\+0\\+","+",formula_levels)
formula_labels <- gsub("~0X\\+","~",formula_labels)
formula_labels <- gsub("~1X\\+","~X+",formula_labels)

cbind(formula_levels,formula_labels)

rbind(varying_intercept[,c("type","formula","x","y")],
      varying_slope[,c("type","formula","x","y")],
      varying_spread[,c("type","formula","x","y")]
      ) %>% mutate(formula=factor(formula,levels=formula_levels,labels=formula_labels)) -> synthetic_data

usethis::use_data(synthetic_data, overwrite = TRUE)

# ggplot(synthetic_data,aes(x=x,y=y,color=formula)) +
#   geom_point(alpha=0.4) +
#   geom_smooth() +
#   facet_wrap(~type) +
#   theme_bw()
#
#
# var_intercept_pwc <- gpct::pairwise_numeric("y",attributes = c("x","formula"),
#                                       data = synthetic_data[which(synthetic_data$type=="varying intercept"),])
# var_slope_pwc <- gpct::pairwise_numeric("y",attributes = c("x","formula"),
#                                   data = synthetic_data[which(synthetic_data$type=="varying slope"),])
# var_spread_pwc <- gpct::pairwise_numeric("y",attributes = c("x","formula"),
#                                         data = synthetic_data[which(synthetic_data$type=="varying spread"),])
#
# var_intercept_results <- gpct::gpct(var_intercept_pwc,explanatory_var = "x",stratum = "formula")
# var_slope_results <- gpct::gpct(var_slope_pwc,explanatory_var = "x",stratum = "formula")
# var_spread_results <- gpct::gpct(var_spread_pwc,explanatory_var = "x",stratum = "formula")
#
# lapply(list(var_intercept=var_intercept_results,
#             var_slope=var_slope_results,
#             var_spread=var_spread_results),function(results){
#
#   cat("\n")
#
#   do.call("rbind",lapply(names(results$strata),function(i){
#     data.frame(strata=i,logOR=results$strata[[i]]$logOR, SE=results$strata[[i]]$logSE)
#   })) %>%
#     mutate(cil=logOR + SE*qnorm(0.05/2),
#            ciu=logOR + SE*qnorm(1-0.05/2),
#     ) %>% print
#
#   print(do.call("c",results$pooled_statistics))
#
# }) %>% invisible()

