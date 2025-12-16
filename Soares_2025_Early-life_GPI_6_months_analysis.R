rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

#======================================================================#
#                                                                      #
#          Costs of early-life GIP infections in spotted hyenas        #
#    repetition of the analysis with a strict criteria (< 6 months)    #
#======================================================================#

##### Packages ####
library(here)
here::here()
library(readr)
library(dplyr)
library(tidyr)
library(meantables)
library(lubridate)
library(lmtest)
library(survival)
library(survminer)
library(ggplot2)
library(ggforestplot)
library(gridExtra)
library(ggiraphExtra)
library(patchwork)
require(ggeffects)
library(broom)
library(scales)

##### Dataset importation and formatting #####

F_early_long  <- read_csv(here("F_early_long.csv"))

# Removal of outlier (fGCM extreme value):

F_early_long <- subset(F_early_long, `f-GCM` < 6000 )

# Repetition of the analysis with a strict criteria (<6months)

F_early_long_6m <- subset(F_early_long,F_early_long$sample_days<= 183)
dim(F_early_long_6m) #51obs
length(unique(F_early_long_6m$ID)) # 51 female individuals
sum(duplicated(F_early_long_6m)) # 0  duplicates

# Normalization of variables to mean 0 and sd 1
F_early_long_6m$`maternal rank` <- scale(F_early_long_6m$`maternal rank`)
F_early_long_6m$`Ancylostoma egg load` <- scale(F_early_long_6m$`Ancylostoma egg load`)
F_early_long_6m$Polyparasitism <- scale(F_early_long_6m$Polyparasitism)
F_early_long_6m$`f-IgA` <- scale(F_early_long_6m$`f-IgA`)
F_early_long_6m$`f-mucin` <- scale(F_early_long_6m$`f-mucin`)
F_early_long_6m$`f-GCM` <- scale(F_early_long_6m$`f-GCM`)

model_survival_6m <- glm (Survival_Ad ~ `maternal rank` +
                            `Ancylostoma egg load` + Polyparasitism +`f-IgA` + `f-mucin` + `f-GCM`,
                          family = "binomial", data = F_early_long_6m)
summary(model_survival_6m)

exp(cbind(Odds_Ratio = coef(model_survival), confint(model_survival)))

# LRT
update_nested <- function(object, formula., ..., evaluate = TRUE){
  update(object = object, formula. = formula., data = object$model, ..., evaluate = evaluate)
}

m1 =  glm (Survival_Ad ~ `maternal rank` +
             `f-IgA` + `f-mucin` + Polyparasitism +  `Ancylostoma egg load` +
             `f-GCM`,
           family = "binomial", data = F_early_long_6m)
m2 = update_nested(m1, . ~ . - `maternal rank`)
m3 = update_nested(m1, . ~ . - `f-IgA`)
m4 = update_nested(m1, . ~ . - `Ancylostoma egg load`)
m5 = update_nested(m1, . ~ . - `maternal rank` - `f-IgA` - `Ancylostoma egg load` - `f-mucin` - Polyparasitism - `f-GCM`)
m6 = update_nested(m1, . ~ . - `f-mucin`)
m7 = update_nested(m1, . ~ . - Polyparasitism)
m8 = update_nested(m1, . ~ . - `f-GCM`)

lrtest(m4,m1) # Ancylostoma egg load
lrtest(m3,m1) # f-IgA
lrtest(m2,m1) # maternal rank
lrtest(m6,m1) # f-mucin
lrtest(m7,m1) # Polyparasitism
lrtest(m8,m1) # f-GCM
