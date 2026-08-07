# siloes.R
# Import siloes data into R
# Reproduce the poisson and mepoisson analyses
# Then run metaGLMM
# May be a bad example because tau-hat=0
# IW 7aug2026

# Set up
library(haven)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(metaGLMM)
setwd("C:/ian/git/metafish/brcancer")
siloes<-read_dta("siloes25_12_meta_pseudo.dta")
siloes$Events = 10 * siloes$Events # makes est tausq=0.3
head(siloes)

# Standard GLM analysis
model <- glm(Events ~ Treat + factor(Study) + offset(log(Pyears)),
             family = poisson, data = siloes)
# coeff & SE of Treat
cbind(
  Estimate = coef(model),
  SE = sqrt(diag(vcov(model)))
)["Treat", ]
#    Estimate          SE 
# -0.31893497  0.04020933 

# Standard GLMM analysis
library(lme4)
model <- glmer(Events ~ Treat + factor(Study) + offset(log(Pyears)) + 
             (0 + Treat | Study),
             family = poisson, data = siloes)
# coeff & SE of Treat
c(fixef(model)["Treat"], sqrt(vcov(model)["Treat", "Treat"]))
#      Treat            
# -0.3680092  0.1181917 
# estimated variance
as.data.frame(VarCorr(model))
#     grp  var1 var2      vcov     sdcor
# 1 Study Treat <NA> 0.3015113 0.5491005


# metaGLMM analysis
# following https://github.com/keisuke-hanada/metaGLMM/readme.md
# Transforming outcomes
dat2 <- data.frame(
  yk = siloes$Events,
  tk = siloes$Pyears,
  Treat = siloes$Treat,
  Study = siloes$Study
) %>%
  mutate(
    thetak = log((pmax(0.5, yk))/tk),
    yk = yk / tk,
    vk = 1 / pmax(0.5,siloes$Events)
  )
# Apply meta-analysis for GLMM
formula <- yk ~ 1 + Treat
dat3 <- dat2 %>% model.frame(formula=formula)
ma.grma1 <- metaGLMM(formula, data=dat3, vi=dat2$vk, ni=dat2$tk, tau2=NA, 
                family=poisson(link="log"), tau2_var=TRUE, ghq_Q=1000)
ma.grma2 <- metaGLMM(formula, data=dat3, vi=dat2$vk, ni=dat2$tk, tau2=NA, 
                family=poisson(link="log"), tau2_var=TRUE, ghq_Q=1000)
summary(ma.grma1)
# Coefficients:
#              Estimate Std. Error z value     Pr(z)    
# (Intercept)  3.851857   0.103852 37.0900 < 2.2e-16 ***
# Treat       -0.340581   0.149440 -2.2790   0.02266 *  
# tau2         0.255674   0.056484  4.5264 5.999e-06 ***
summary(ma.grma2)
# Coefficients:
#              Estimate Std. Error z value     Pr(z)    
# (Intercept)  3.851867   0.103851 37.0903 < 2.2e-16 ***
# Treat       -0.340589   0.149439 -2.2791   0.02266 *  
# tau2         0.255670   0.056483  4.5265 5.997e-06 ***

### confidence interval by profile likelihood
ci.pl <- confint_PL(ma.grma2)
ci.pl
#                  lower       upper
# (Intercept)  3.6436247  4.05941685
# Treat       -0.6413195 -0.04373964

### confidence interval by profile likelihood with simple Bartlett correction
ci.plsbc <- confint_SBC(ma.grma2)
ci.plsbc
#                  lower       upper
# (Intercept)  3.6403621  4.06267226
# Treat       -0.6460376 -0.03913956
