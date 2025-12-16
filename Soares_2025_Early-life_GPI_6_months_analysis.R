<<<<<<< HEAD
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

###### Analyses ######

###### 1.1. Survival to adulthood ######
# Repetition of the analysis with a strict criteria (<6months)

# Sample size:
F_early_long_6m <- subset(F_early_long,F_early_long$sample_days<= 183)
dim(F_early_long_6m) #51obs
length(unique(F_early_long_6m$ID)) # 51 female individuals
sum(duplicated(F_early_long_6m)) # 0  duplicates

# Checking formatting of response variable
is.factor(F_early_long_6m$Survival_Ad) #FALSE
F_early_long_6m$Survival_Ad <- as.factor(F_early_long_6m$Survival_Ad)
contrasts(F_early_long_6m$Survival_Ad)

F_early_long_6m$Survival_Ad <- droplevels(F_early_long_6m$Survival_Ad)
levels(F_early_long_6m$Survival_Ad)

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

##### 1.2. Age at First Reproduction (AFR) #####

# Removal individuals that didn´t survived until adulthood
# Restrict analysis to 6months of age mark
F_early_long_AFR <- subset(F_early_long, Survival_Ad == "Yes" )
F_early_long_AFR_6m <- subset(F_early_long_AFR,F_early_long_AFR$sample_days<= 183)

# Create Event variable:  reproduced or not
F_early_long_AFR_6m <- F_early_long_AFR_6m %>%
  mutate(Reproduced = case_when(AFR >= 0  ~ 1,
                                is.na(AFR) ~ 0))

# Create variable for Time to Event (YEAR)
#females with AFR -> AFR
#females never reproduced and dead -> age of death
#female never reproduced still alive -> current date of dataset - birthdate
current_date = as.Date("2023-03-23", format = "%Y-%m-%d")
library(lubridate)
int <- interval(ymd(F_early_long_AFR_6m$birthdate), ymd(F_early_long_AFR_6m$deathdate))
F_early_long_AFR_6m$age_death <- time_length(int, "year")

int1 <- interval(ymd(F_early_long_AFR_6m$birthdate), ymd(current_date))
F_early_long_AFR_6m$current_age <- time_length(int1, "year")

F_early_long_AFR_6m<- F_early_long_AFR_6m%>%
  mutate(tafr = case_when(AFR >=0 ~ AFR,
                          is.na(AFR) & age_death >= 0 ~ age_death,
                          is.na(AFR) & is.na(age_death) ~ current_age))

# Cox proportional hazards model:

# # Normalization of variables to mean 0 and sd 1
F_early_long_AFR_6m$`Ancylostoma egg load` <- scale(F_early_long_AFR_6m$`Ancylostoma egg load`)
F_early_long_AFR_6m$`f-IgA` <- scale(F_early_long_AFR_6m$`f-IgA`)
F_early_long_AFR_6m$`f-GCM` <- scale(F_early_long_AFR_6m$`f-GCM`)
F_early_long_AFR_6m$`maternal rank 6m` <- scale(F_early_long_AFR_6m$`maternal rank 6m`)

AFR_object_6m <- Surv(time = F_early_long_AFR_6m$tafr, event = F_early_long_AFR_6m$Reproduced)
AFR_object_6m
AFR_model_6m <- coxph(AFR_object_6m ~ `maternal rank 6m` +  `Ancylostoma egg load` +  `f-IgA` + `f-GCM` , data = F_early_long_AFR_6m)
summary(AFR_model_6m)

# LRT from anova
m1 <- coxph(AFR_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` + `f-IgA` + `f-GCM` , data = F_early_long_AFR_6m)
m0 <- update(m1, . ~ 1)
m2 <- coxph(AFR_object_6m ~  `Ancylostoma egg load` + `f-IgA` + `f-GCM`, data =F_early_long_AFR_6m)
m3 <- coxph(AFR_object_6m ~  `maternal rank 6m` +  `f-IgA` +  `f-GCM`, data = F_early_long_AFR_6m)
m4 <- coxph(AFR_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` +  `f-GCM` , data = F_early_long_AFR_6m)
m5 <- coxph(AFR_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` +  `f-IgA` , data = F_early_long_AFR_6m)

anova(m2,m1) #maternal rank
anova(m3,m1) #Ancylostoma egg load
anova(m4,m1) #f-IgA
anova(m5,m1) #f-GCM

###### 1.3. Longevity ######
=======
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

###### Analyses ######

###### 1.1. Survival to adulthood ######
# Repetition of the analysis with a strict criteria (<6months)

# Sample size:
F_early_long_6m <- subset(F_early_long,F_early_long$sample_days<= 183)
dim(F_early_long_6m) #51obs
length(unique(F_early_long_6m$ID)) # 51 female individuals
sum(duplicated(F_early_long_6m)) # 0  duplicates

# Checking formatting of response variable
is.factor(F_early_long_6m$Survival_Ad) #FALSE
F_early_long_6m$Survival_Ad <- as.factor(F_early_long_6m$Survival_Ad)
contrasts(F_early_long_6m$Survival_Ad)

F_early_long_6m$Survival_Ad <- droplevels(F_early_long_6m$Survival_Ad)
levels(F_early_long_6m$Survival_Ad)

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

##### 1.2. Age at First Reproduction (AFR) #####

# Removal individuals that didn´t survived until adulthood
# Restrict analysis to 6months of age mark
F_early_long_AFR <- subset(F_early_long, Survival_Ad == "Yes" )
F_early_long_AFR_6m <- subset(F_early_long_AFR,F_early_long_AFR$sample_days<= 183)

# Create Event variable:  reproduced or not
F_early_long_AFR_6m <- F_early_long_AFR_6m %>%
  mutate(Reproduced = case_when(AFR >= 0  ~ 1,
                                is.na(AFR) ~ 0))

# Create variable for Time to Event (YEAR)
#females with AFR -> AFR
#females never reproduced and dead -> age of death
#female never reproduced still alive -> current date of dataset - birthdate
current_date = as.Date("2023-03-23", format = "%Y-%m-%d")
library(lubridate)
int <- interval(ymd(F_early_long_AFR_6m$birthdate), ymd(F_early_long_AFR_6m$deathdate))
F_early_long_AFR_6m$age_death <- time_length(int, "year")

int1 <- interval(ymd(F_early_long_AFR_6m$birthdate), ymd(current_date))
F_early_long_AFR_6m$current_age <- time_length(int1, "year")

F_early_long_AFR_6m<- F_early_long_AFR_6m%>%
  mutate(tafr = case_when(AFR >=0 ~ AFR,
                          is.na(AFR) & age_death >= 0 ~ age_death,
                          is.na(AFR) & is.na(age_death) ~ current_age))

# Cox proportional hazards model:

# # Normalization of variables to mean 0 and sd 1
F_early_long_AFR_6m$`Ancylostoma egg load` <- scale(F_early_long_AFR_6m$`Ancylostoma egg load`)
F_early_long_AFR_6m$`f-IgA` <- scale(F_early_long_AFR_6m$`f-IgA`)
F_early_long_AFR_6m$`f-GCM` <- scale(F_early_long_AFR_6m$`f-GCM`)
F_early_long_AFR_6m$`maternal rank 6m` <- scale(F_early_long_AFR_6m$`maternal rank 6m`)

AFR_object_6m <- Surv(time = F_early_long_AFR_6m$tafr, event = F_early_long_AFR_6m$Reproduced)
AFR_object_6m
AFR_model_6m <- coxph(AFR_object_6m ~ `maternal rank 6m` +  `Ancylostoma egg load` +  `f-IgA` + `f-GCM` , data = F_early_long_AFR_6m)
summary(AFR_model_6m)

# LRT from anova
m1 <- coxph(AFR_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` + `f-IgA` + `f-GCM` , data = F_early_long_AFR_6m)
m0 <- update(m1, . ~ 1)
m2 <- coxph(AFR_object_6m ~  `Ancylostoma egg load` + `f-IgA` + `f-GCM`, data =F_early_long_AFR_6m)
m3 <- coxph(AFR_object_6m ~  `maternal rank 6m` +  `f-IgA` +  `f-GCM`, data = F_early_long_AFR_6m)
m4 <- coxph(AFR_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` +  `f-GCM` , data = F_early_long_AFR_6m)
m5 <- coxph(AFR_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` +  `f-IgA` , data = F_early_long_AFR_6m)

anova(m2,m1) #maternal rank
anova(m3,m1) #Ancylostoma egg load
anova(m4,m1) #f-IgA
anova(m5,m1) #f-GCM

###### 1.3. Longevity ######
>>>>>>> d4726b1cba48d50c732deeb512145d57235518cc
