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
library(patchwork)
library(glmmTMB)  # for poisson regression.

source("./scripts/theme-opts.R")

###############################################################################
## Experiment 2 (Fig 1 c and d)
###############################################################################


## Raw  data survivorship numbers by vial and by day:
cd <- read_csv("./data/fig1cd-survival.csv")

## table assigning vials to treatments/source population
vials <- read_csv("./data/fig1cd-vials.csv") %>% select(-graph)
vials$treatment <- factor(vials$treatment, levels = c("Females Only",
                                                      "X Spermless males",
                                                      "X Oregon R males"))

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


###############################################################################
## Experiment 1 (Fig 1 a and b)
###############################################################################

## Raw  data survivorship numbers by vial and by day:
b <- read_csv("./data/fig1b-survival.csv")

## table assigning vials to treatments/source population
vials_ab <- read_csv("./data/fig1ab-vials.csv") %>% select(-source_pop)

## Get raw survival data in long format with n_alive as explicit column:
b_raw <- pivot_longer(b, -Week, names_to="vial", values_to="n_alive")
nrow(b_raw)
# 160

## We need to turn the data showing number alive into a data frame showing how
## many died each time step. I'll subtract using lag within group. All vials have
# 50 flies at day 0, so I can hard code that.
b_surv <- b_raw %>% arrange(vial, Week) %>% group_by(vial) %>%
  mutate(mort=lag(n_alive, default=50) - n_alive)

# Now we need to turn this into one row per fly as is expected by the survival
# package, that means "uncounting" the data summarized to vial. 
b_surv_long <- b_surv %>% select(vial, Week, mort) %>%
  uncount(weights=mort)
nrow(b_surv_long)


b_surv_long$status <- 1 # all our events are death


## 2025-01-24: not needed, Lewis had data for week 9 when all flies were deda.
## added now.
## Hacky but I will add in living flies at end of experiment by
## hand here: living <- b_raw %>% filter(Week==8 & n_alive > 0) %>%
## uncount(weight=n_alive) living$status <- 0 # means alive b_surv_long <-
## bind_rows(b_surv_long, living)

## Now merge this with the population and treatment data associated with each
## vial:
b_surv_final <- left_join(b_surv_long, vials_ab)
nrow(b_surv_long)
# 1000


## Now our data is ready for the survival package. We read it in as a "Surv" object
b_surv_obj <- Surv(time=b_surv_final$Week, b_surv_final$status)

# Let's make sure our Kaplan-Meier plots look like the ms figure:
KM_fit_b <- survfit(b_surv_obj ~ treatment, data=b_surv_final)
KM_fit_b
#plot(KM_fit, xlab = "Time to Death (days)", ylab = "Survival Probability")
ggsurvplot(KM_fit_b)
## OK, that looks like our six curves
ggsave("./results/KM_plot_b.pdf")

ggsurvplot( survfit(cd_surv_obj ~ vial, data=cd_surv_final))

## Now for the mixed effect Cox model
coxmod_b <- coxme(b_surv_obj ~ treatment + (1 | vial), data=b_surv_final)
coxmod_b
## Cox mixed-effects model fit by maximum likelihood
##   Data: b_surv_final
##   events, n = 995, 1000
##   Iterations= 7 38 
##                     NULL Integrated    Fitted
## Log-likelihood -5907.341  -5531.546 -5513.627

##                    Chisq    df p    AIC    BIC
## Integrated loglik 751.59  2.00 0 747.59 737.78
##  Penalized loglik 787.43 13.71 0 760.01 692.78

## Model:  b_surv_obj ~ treatment + (1 | vial) 
## Fixed coefficients
##                                coef exp(coef)  se(coef)     z p
## treatmentX Spermless Males 2.664505  14.36083 0.1492842 17.85 0

## Random effects
##  Group Variable  Std Dev    Variance  
## vial  Intercept 0.22096480 0.04882544

## Anova table:
car::Anova(coxmod_b)

## Analysis of Deviance Table (Type II tests)

## Response: b_surv_obj
##           Df  Chisq Pr(>Chisq)    
## treatment  1 318.57  < 2.2e-16 ***
## ---
##   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1






###############################################################################
## Figure 1 a,b,c,d
###############################################################################

## Figure theme elements:
color_scale <- scale_color_brewer(palette = "Dark2")
error_width <- 0.3
line_width <- 0.55


# Fig 1 c and d:
surv_fig_data <- left_join(cd_raw, vials)  %>% mutate(survival = n_alive) #/50

mean_fig_data <- surv_fig_data %>%
  group_by(Day, source_pop, treatment) %>%
  summarize(m_survival = mean(survival), se = sd(survival)/sqrt(sum(!is.na(survival))),
            lower = max(m_survival - se, 0), upper = m_survival + se)


## Fig 1c:
fig1c <- ggplot(filter(mean_fig_data, source_pop=="Parthenogenic females"),
                aes(x=Day, y=m_survival, color=treatment)) +
  ggtitle("c) Parthenogenic females") +
  color_scale +
  geom_point(size=0.5) +
  geom_line(linewidth=line_width) +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=0,
                linewidth=error_width,
                alpha=0.7) +
  pubtheme.nogridlines +
  xlab("Days after mating") +
  ylab("# of females alive per vial") +
  ylim(c(0,50)) +
  theme(legend.position=c(0.7, 0.85),
        legend.title=element_blank(),
        legend.text = element_text(family=fontfamily, size=smsize-2))
ggsave("./results/fig1c.pdf", plot=fig1c, width=col1, height=col1, units="cm")

fig1d <- ggplot(filter(mean_fig_data, source_pop=="Oregon R females"),
                aes(x=Day, y=m_survival, color=treatment)) +
  ggtitle("d) Oregon R females") +
  color_scale +
  geom_point(size=0.5) +
  geom_line(linewidth=line_width) +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=0,
                linewidth=error_width,
                alpha=0.7) +
  pubtheme.nogridlines +
  xlab("Days after mating") + ylab("# of females alive per vial") +
  ylim(c(0,50)) +
  theme(legend.position="none")
  ## theme(legend.position=c(0.7, 0.85),
  ##       legend.title=element_blank(),
  ##       legend.text = element_text(family=fontfamily, size=smsize-2))
ggsave("./results/fig1d.pdf", plot=fig1d, width=col1, height=col1, units="cm")


### Fig 1b 
surv_fig_data_b <- left_join(b_raw, vials_ab)  %>% mutate(survival = n_alive) #/50

mean_fig_data_b <- surv_fig_data_b %>%
  group_by(Week, treatment) %>%
  summarize(m_survival = mean(survival), se = sd(survival)/sqrt(sum(!is.na(survival))),
            lower = max(m_survival - se, 0), upper = m_survival + se)

fig1b <- ggplot(mean_fig_data_b, aes(x=Week, y=m_survival, color=treatment)) +
  xlim(c(0.5,8)) +
  ggtitle("b) Parthenogenic females") +
  color_scale +
  geom_point(size=0.5) +
  geom_line(linewidth=line_width) +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=0,
                linewidth=error_width,
                alpha=0.7) +
  pubtheme.nogridlines +
  xlab("Weeks after mating") + ylab("# of females alive per vial") +
  ylim(c(0,50)) +
   theme(legend.position="none")
  ##theme(legend.position=c(0.72, 0.88),
  ##       legend.title=element_blank(),
  ##       legend.text = element_text(family=fontfamily, size=smsize-2))
ggsave("./results/fig1b.pdf", plot=fig1b, width=col1, height=col1, units="cm")


## Fig 1a

fecund <- read_csv("./data/fig1a-fecund.csv")
fecund <- pivot_longer(fecund, -Week, names_to="vial", values_to="n_f1")
nrow(fecund) # 160
fecund <- left_join(fecund, vials_ab)

mean_fecund <- fecund %>%
  group_by(Week, treatment) %>%
  summarize(mean_f1 = mean(n_f1), se = sd(n_f1)/sqrt(sum(!is.na(n_f1))),
            lower = max(mean_f1 - se, 0), upper = mean_f1 + se)


fig1a <- ggplot(mean_fecund, aes(Week, mean_f1, color=treatment)) +
   xlim(c(0.5,8)) +
  ggtitle("a) F1 of Parthenogenic females") +
  color_scale +
  geom_point(size=0.5) +
  geom_line(linewidth=line_width) +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=0,
                linewidth=error_width,
                alpha=1) +
  pubtheme.nogridlines +
  xlab("Weeks after mating") +
  ylab("# of F1 adults / 50 females") +
  theme(legend.position=c(0.72, 0.88),
        legend.title=element_blank(),
        legend.text = element_text(family=fontfamily, size=smsize-2))
ggsave("./results/fig1a.pdf", plot=fig1a, width=col1, height=col1, units="cm")

## Compare mean times to F1 adult in weeks. Observations are F1s that reached adult.
fecund_ind <- uncount(fecund, weights=n_f1)
fecund_ind
fecund_mod <- glmmTMB(Week ~ treatment + (1 | vial), data = fecund_ind, family = poisson)
summary(fecund_mod)
##  Family: poisson  ( log )
## Formula:          Week ~ treatment + (1 | vial)
## Data: fecund_ind

##      AIC      BIC   logLik deviance df.resid 
##    112.3    116.8    -53.1    106.3       30 

## Random effects:

## Conditional model:
##  Groups Name        Variance  Std.Dev.
##  vial   (Intercept) 2.015e-10 1.42e-05
## Number of obs: 33, groups:  vial, 18

## Conditional model:
##                            Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                  1.3041     0.1195  10.911  < 2e-16 ***
## treatmentX Spermless Males  -0.6109     0.2236  -2.732  0.00629 ** 
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
car::Anova(fecund_mod)

## Analysis of Deviance Table (Type II Wald chisquare tests)

## Response: Week
##            Chisq Df Pr(>Chisq)   
## treatment 7.4642  1   0.006294 **
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## > fecund_ind



## s_fecund <- fecund %>% group_by(vial, treatment) %>% summarize(total_f1 = sum(n_f1))
## lm(total_f1 ~  treatment, data=s_fecund)
## summary(fecund_sum_mod)
mean(fecund_ind$Week[fecund_ind$treatment=="Females Only"]) # 3.7
mean(fecund_ind$Week[fecund_ind$treatment=="X Spermless Males"]) # 2.0


# Total # f1 adults per vial produced not significantly different across
# treatments (p = 0.2691
s_fecund <- fecund %>% group_by(vial, treatment) %>% summarize(total_f1 = sum(n_f1))
fecund_sum_mod <- lm(total_f1 ~  treatment, data=s_fecund)
summary(fecund_sum_mod)





#### combined figure 1
fig1 <- fig1a + fig1b + fig1c + fig1d +
  plot_layout(nrow = 2, byrow = FALSE)
#  plot_annotation(tag_levels = 'A') &
#  theme(plot.margin = unit(c(4, 4, 4, 4), "pt"),
#        plot.tag.position = c(0, 1),
#        plot.tag = element_text(hjust = -0.5, vjust = 0.3))
fig1

ggsave("./results/fig1.pdf", plot=fig1,
       width=col2, height=col2, unit="cm")
