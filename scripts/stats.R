## Dylan Schwilk
## 2025-01-23

## Script to compare survivorship curves across two drosphila populations and
## three treatments.

## packages needed
library(readr)
library(tidyr)
library(dplyr)
library(survival)
library(coxme) # for mixed effect Cox model
library(survminer) # for prettier plots using ggplot2
# library(car)

source("./scripts/theme-opts.R")

## Raw  data survivorship numbers by vial and by day:
cd <- read_csv("./data/fig1cd-survival.csv")

## table assigning vials to treatments/source population
vials <- read_csv("./data/fig1cd-vials.csv") %>% select(-graph)

## Get raw survival data in long format with n_alive as explicit column:
cd_raw <- pivot_longer(cd, -Day, names_to="vial", values_to="n_alive")
nrow(cd_raw)
# 888 = 24 vials by 37 times, good

## We need to turn the data showing number alive into a data frame showing how
## many died each time step. I'll subtract using lag within group. All vials have
# 50 flies at day 0, so I can hard code that.
cd_surv <- cd_raw %>% arrange(vial, Day) %>% group_by(vial) %>%
  mutate(mort=lag(n_alive, default=50)-n_alive)

# Now we need to turn this into one row per fly as is expected by the survival
# package, that means "uncounting" the data summarized to vial. 
cd_surv_long <- cd_surv %>% select(vial, Day, mort) %>%
  uncount(weights=mort)
nrow(cd_surv_long)
# 1200 = 24 vials * 50 flies, good.

## Now merge this with the population and treatment data associated with each
## vial:
cd_surv_final <- left_join(cd_surv_long, vials)
cd_surv_final$status <- 1 # all our events are death

## Now our data is ready for the survival package. We read it in as a "Surv" object
cd_surv_obj <- Surv(time=cd_surv_final$Day, cd_surv_final$status)

# Let's make sure our Kaplan-Meier plots look like the ms figure:
KM_fit <- survfit(cd_surv_obj ~ source_pop + treatment, data=cd_surv_final)
KM_fit
#plot(KM_fit, xlab = "Time to Death (days)", ylab = "Survival Probability")
ggsurvplot(KM_fit)
## OK, that looks like our six curves
ggsave("./results/KM_plot_all.pdf")

ggsurvplot( survfit(cd_surv_obj ~ vial, data=cd_surv_final))

## Now for the mixed effect Cox model
coxmod <- coxme(cd_surv_obj ~ source_pop * treatment + (1 | vial), data=cd_surv_final)
coxmod
## Cox mixed-effects model fit by maximum likelihood
##   Data: cd_surv_final
##   events, n = 1200, 1200
##   Iterations= 11 71 
##                     NULL Integrated   Fitted
## Log-likelihood -7312.556  -6738.614 -6687.01

##                     Chisq    df p     AIC     BIC
## Integrated loglik 1147.88  6.00 0 1135.88 1105.34
##  Penalized loglik 1251.09 22.47 0 1206.14 1091.75

## Model:  cd_surv_obj ~ source_pop * treatment + (1 | vial) 
## Fixed coefficients
##                                                                 coef   exp(coef)  se(coef)     z       p
## source_popParthenogenic females                             2.150383  8.58814400 0.6057433  3.55 3.9e-04
## treatmentX Oregon R males                                   2.358558 10.57568822 0.6076604  3.88 1.0e-04
## treatmentX Spermless males                                  4.083867 59.37461859 0.6106077  6.69 2.3e-11
## source_popParthenogenic females:treatmentX Oregon R males  -2.244307  0.10600098 0.8542735 -2.63 8.6e-03
## source_popParthenogenic females:treatmentX Spermless males -2.750554  0.06389245 0.8548252 -3.22 1.3e-03

## Random effects
##  Group Variable  Std Dev   Variance 
## vial  Intercept 0.8395459 0.7048374


## Anova table:
car::Anova(coxmod)
## Analysis of Deviance Table (Type II tests)

## Response: cd_surv_obj
##                      Df   Chisq Pr(>Chisq)    
## source_pop            1  1.8969   0.168424    
## treatment             2 39.0692  3.283e-09 ***
## source_pop:treatment  2 11.7117   0.002863 ** 
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



surv_fig_data <- left_join(cd_raw, vials) %>% mutate(survival = n_alive/50)


mean_fig_data <- surv_fig_data %>%
  group_by(Day, source_pop, treatment) %>%
  summarize(survival = mean(survival), sd=sd(survival))


stat_sum_single <- function(fun, geom="point", ...) {
  stat_summary(fun.y=fun, geom=geom, size = 3, ...)
}

ggplot(mean_fig_data, aes(x=Day, y=survival, color=treatment)) +
#  geom_point(size=0.5) +
  geom_line(linewidth=0.6) +
  stat_summary(fun.data="mean_sdl",
               fun.args = list(mult = 1),
               position = position_jitter(width=0.2, height=0),
               size=0.1,
               linewidth=0.6,
               alpha=0.8,
               data=surv_fig_data) +
#geom_errorbar(aes(ymin=survival - sd, ymax=survival + sd), width=.8) +
  facet_grid(source_pop ~ .) +
  pubtheme.nogridlines +
  xlab("Time (d)") + ylab("Surviving proportion") +
  theme(legend.position=c(0.75, 0.45),
        legend.title=element_blank(),
        legend.text = element_text(family=fontfamily, size=smsize-2))
#       legend.key.height=unit(smsize,"pt"))
  

ggsave("./results/fig1cd.pdf", width=col1, height=col2, units="cm")
