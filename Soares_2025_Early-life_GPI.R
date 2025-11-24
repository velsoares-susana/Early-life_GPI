rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

#############################################################
####       Early adversity - long term effects         #####
############################################################
#packages#
install.packages("here")
install.packages("tidyr")
install.packages("meantables")
install.packages("lubridate")
install.packages("lmtest")
install.packages("survival")
install.packages("survminer")
install.packages("ggforestplot")
install.packages("broom")


# Load required packages
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


#if need be
ls() # # to see if there are any objects present loaded in the working environment
rm(list=ls()) # # clearing the environment before running the script

#database to work on
#F_early_long

dim(F_early_long) #69 obs
#removal of cortisol outlier
F_early_long <- subset(F_early_long, `f-GCM` < 6000 )

#descriptives

F_early_long %>%
  mean_table(AFR)
F_early_long %>%
  mean_table(Longevity)
F_early_long %>%
  mean_table(sample_years)

####plot sampling/lab effort -  Fig.S1 ####

#samples completed for all assays
F_early_long$all_sample <- case_when(F_early_long$`f-mucin` >= 0 & F_early_long$Polyparasitism >= 0 & F_early_long$`f-IgA` >=0 & F_early_long$`f-GCM` >=0 ~ F_early_long$sample_days)
#delimiting lifelines
F_early_long$age_death<- difftime(F_early_long$deathdate ,F_early_long$birthdate , units = c("days"))
d_numeric <- sapply(F_early_long, class) == "difftime"
F_early_long[d_numeric] <- lapply(F_early_long[d_numeric], as.numeric)

F_early_long$age_death_c <- case_when(F_early_long$age_death > 365  ~ 400,
                                      F_early_long$age_death < 365  ~ F_early_long$age_death,
                                      TRUE ~ 400)

sampling_distribution_all <- ggplot(F_early_long, aes(ID, sample_days)) +
  geom_segment((aes(x=ID, xend = ID, y= 0, yend = sample_days)), colour ="lightgrey") +
  geom_segment((aes(x=ID, xend = ID, y= sample_days, yend = age_death_c)), colour ="lightgrey")+
  geom_point(aes(x = ID, y = all_sample), colour = "black")+
  geom_hline(yintercept=183, linetype="dashed") +
  labs (x = "Females sampled (2010-2017)", y= "Age at sampling (days)") +
  scale_y_continuous(expand = c(0, 0),  breaks = c(0, 183, 365),labels = c(0,"6m", "1y"), limits = c(0, 400)) +
  coord_flip() +
  theme_classic()
sampling_distribution_final_all <- sampling_distribution_all +
  theme (axis.text.y = element_blank(), axis.ticks.y = element_blank())
sampling_distribution_final_all





##################################################################################################
###        Models
##############################################
################################################################################################################
##  1.1 Survival to adulthood  ##
#################################################################################################################
#sample size#

dim(F_early_long) #68 obs
ggplot(F_early_long, aes (x = Survival_Ad)) +
  geom_bar() +
  labs (x = "Survival_Ad", y = "count")
table(F_early_long$Survival_Ad) # NO 21 Yes 48

#assumptions
# checking formatting of variable
is.factor(F_early_long$Survival_Ad) #FALSE
F_early_long$Survival_Ad <- as.factor(F_early_long$Survival_Ad)
contrasts(F_early_long$Survival_Ad)

F_early_long$Survival_Ad <- droplevels(F_early_long$Survival_Ad)
levels(F_early_long$Survival_Ad)

##data distribution plots - Fig.S2 ##
#independent variables
# maternal rank
plot_rank <- ggplot(F_early_long, aes(`maternal rank`)) +
  geom_histogram(color = "#000000", fill = "grey") +
  labs(,
       x = "maternal social rank",
       y = "Count"
  ) +
  theme_classic()
#f-IgA
plot_IgA <- ggplot(F_early_long, aes(`f-IgA`)) +
  geom_histogram(color = "#000000", fill = "grey") +
  labs(,
       x = "f-IgA levels",
       y = "Count"
  ) +
  theme_classic()
#f-mucin
plot_mucin<- ggplot(F_early_long, aes(`f-mucin`)) +
  geom_histogram(color = "#000000", fill = "grey") +
  labs(,
       x = "f-mucin levels",
       y = "Count"
  ) +
  theme_classic()
#Ancylostoma egg load
plot_Ancy <- ggplot(F_early_long, aes(`Ancylostoma egg load`)) +
  geom_histogram(color = "#000000", fill = "grey") +
  labs(,
       x = "Ancylostoma egg load",
       y = "Count"
  ) +
  theme_classic()
#Polyparasitism
plot_rich <- ggplot(F_early_long, aes(Polyparasitism)) +
  geom_histogram(color = "#000000", fill = "grey") +
  labs(,
       x = "Polyparasitism without presence of Ancylostoma sp.",
       y = "Count"
  ) +
  theme_classic()
#f-GCM
plot_cort <- ggplot(F_early_long, aes(`f-GCM`)) +
  geom_histogram(color = "#000000", fill = "grey") +
  labs(,
       x = "f-GCM levels",
       y = "Count"
  ) +
  theme_classic()

## agregate plots
Sampling_plot <- grid.arrange (plot_rank, plot_Ancy,plot_rich, plot_IgA, plot_mucin, plot_cort, nrow=3)
###

##########################################################
## glm ##
F_early_long_Surv <- F_early_long
#normalization of variables to mean 0 and sd 1
F_early_long_Surv$`maternal rank` <- scale(F_early_long_Surv$`maternal rank`)
F_early_long_Surv$`Ancylostoma egg load` <- scale(F_early_long_Surv$`Ancylostoma egg load`)
F_early_long_Surv$Polyparasitism <- scale(F_early_long_Surv$Polyparasitism)
F_early_long_Surv$`f-IgA` <- scale(F_early_long_Surv$`f-IgA`)
F_early_long_Surv$`f-mucin` <- scale(F_early_long_Surv$`f-mucin`)
F_early_long_Surv$`f-GCM` <- scale(F_early_long_Surv$`f-GCM`)

#model
model_survival <- glm (Survival_Ad ~ `maternal rank` +
                         `Ancylostoma egg load` + Polyparasitism + `f-IgA` + `f-mucin` + `f-GCM`,
                       family = "binomial", data = F_early_long_Surv)
summary(model_survival)

exp(cbind(Odds_Ratio = coef(model_survival), confint(model_survival)))

#################
####### LRT
update_nested <- function(object, formula., ..., evaluate = TRUE){
  update(object = object, formula. = formula., data = object$model, ..., evaluate = evaluate)
}

m1 =  glm (Survival_Ad ~ `maternal rank` +
             `f-IgA` + `f-mucin` + Polyparasitism +  `Ancylostoma egg load` +
             `f-GCM`,
           family = "binomial", data = F_early_long_Surv)
m2 = update_nested(m1, . ~ . - `maternal rank`)
m3 = update_nested(m1, . ~ . - `f-IgA`)
m4 = update_nested(m1, . ~ . - `Ancylostoma egg load`)
m5 = update_nested(m1, . ~ . - `maternal rank` - `f-IgA` - `Ancylostoma egg load` - `f-mucin` - Polyparasitism - `f-GCM`)
m6 = update_nested(m1, . ~ . - `f-mucin`)
m7 = update_nested(m1, . ~ . - Polyparasitism)
m8 = update_nested(m1, . ~ . - `f-GCM`)


lrtest(m5,m1) #intercept
lrtest(m4,m1) # Ancylostoma egg load
lrtest(m3,m1) #f-IgA
lrtest(m2,m1) # maternal rank
lrtest(m6,m1) # f-mucin
lrtest(m7,m1) # Polyparasitism
lrtest(m8,m1) # f-GCM

##########################################################
##plot of significant effects - Fig. 1
F_early_long_Surv_plot <- F_early_long
is.factor(F_early_long_Surv_plot$Survival_Ad) #FALSE
F_early_long_Surv_plot$Survival_Ad <- as.factor(F_early_long_Surv_plot$Survival_Ad)
contrasts(F_early_long_Surv_plot$Survival_Ad)

F_early_long_Surv_plot$Survival_Ad <- droplevels(F_early_long_Surv_plot$Survival_Ad)
levels(F_early_long_Surv_plot$Survival_Ad)
library(dplyr)
F_early_long_Surv_plot <- F_early_long_Surv_plot %>%
  rename(maternal_rank = `maternal rank`,
         Ancylostoma_egg_load = `Ancylostoma egg load`,
         f_IgA = `f-IgA`, f_mucin = `f-mucin`,
         f_GCM = `f-GCM`)

model_survival <- glm (Survival_Ad ~ maternal_rank +
                         Ancylostoma_egg_load + Polyparasitism + f_IgA + f_mucin + f_GCM,
                       family = "binomial", data = F_early_long_Surv_plot)
summary(model_survival)

# Maternal_plot
Maternal_plot <- ggpredict(model_survival, terms = c("maternal_rank[all]")) |>
  plot(colors = "darkorchid4", show_data = TRUE, jitter = 0.05) +
  labs(subtitle = "a)",
       title = NULL,
       x = "maternal rank",
       y = "Predicted probabilities of Survival to Adulthood") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.subtitle = element_text(face = "bold", size = 18))
Maternal_plot

# IgA_plot
IgA_plot <- ggpredict(model_survival, terms = c("f_IgA[all]")) |>
  plot(colors = "#C71585", show_data = TRUE, jitter = 0.05) +
  labs(subtitle = "b)",
       title = NULL,
       x = "f-IgA level",
       y = NULL) +
  scale_y_continuous(labels = NULL) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.subtitle = element_text(face = "bold", size = 18))

# Ancy_plot
Ancy_plot <- ggpredict(model_survival, terms = c("Ancylostoma_egg_load[all]")) |>
  plot(colors = "#ADD8E6", show_data = TRUE, jitter = 0.05) +
  labs(subtitle = "c)",
       title = NULL,
       x = "Ancylostoma egg load",
       y = NULL) +
  scale_y_continuous(labels = NULL) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.subtitle = element_text(face = "bold", size = 18))

# combine the plots into a single row
Survival_plot <- Maternal_plot + IgA_plot + Ancy_plot
Survival_plot

#save file
ggsave(
  filename = "Survival.png",
  plot = Survival_plot,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)


#################################################################################
##  Confirmation at 6m                                  ##
################################################################################

# restrict dataset to observations within 0.5 years of living
F_early_long_6m <- subset(F_early_long,F_early_long$sample_days<= 183)
dim(F_early_long_6m) #51obs
length(unique(F_early_long_6m$ID)) # 51 female individuals
sum(duplicated(F_early_long_6m)) # 0  duplicates


#################################################
## model ###
################################################
#normalization of variables to mean 0 and sd 1
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

###############
#LRT
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

lrtest(m5,m1)
lrtest(m4,m1) # Ancylostoma egg load
lrtest(m3,m1) # f-IgA
lrtest(m2,m1) # maternal rank
lrtest(m6,m1) # f-mucin
lrtest(m7,m1) # Polyparasitism
lrtest(m8,m1) # f-GCM


##################################################################################
## 1.2. AFR ###
#################################################################################
###########################################
#removed individuals that didn´t survived until adulthood

F_early_long_AFR <- subset(F_early_long, Survival_Ad == "Yes" )


#create event variable:  reproduced or not
F_early_long_AFR <- F_early_long_AFR %>%
  mutate(Reproduced = case_when(AFR >= 0  ~ 1,
                                is.na(AFR) ~ 0))

#create variables for time to event variable (YEAR)
#females with AFR -> AFR
#females never reproduced and dead -> age of death
#female never reproduced still alive -> current date of dataset - birthdate
current_date = as.Date("2023-03-23", format = "%Y-%m-%d")
library(lubridate)
int <- interval(ymd(F_early_long_AFR$birthdate), ymd(F_early_long_AFR$deathdate))
F_early_long_AFR$age_death <- time_length(int, "year")

int1 <- interval(ymd(F_early_long_AFR$birthdate), ymd(current_date))
F_early_long_AFR$current_age <- time_length(int1, "year")


# create time to event variable
F_early_long_AFR<- F_early_long_AFR%>%
  mutate(tafr = case_when(AFR >=0 ~ AFR,
                          is.na(AFR) & age_death >= 0 ~ age_death,
                          is.na(AFR) & is.na(age_death) ~ current_age))


# Fit the Cox proportional hazard model
AFR_object <- Surv(time = F_early_long_AFR$tafr, event = F_early_long_AFR$Reproduced)
AFR_object

# Fit a Cox proportional hazards model with normalized variables
# scale variables to mean 0 and sd 1 #
F_early_long_AFR$`Ancylostoma egg load` <- scale(F_early_long_AFR$`Ancylostoma egg load`)
F_early_long_AFR$`f-IgA` <- scale(F_early_long_AFR$`f-IgA`)
F_early_long_AFR$`f-GCM` <- scale(F_early_long_AFR$`f-GCM`)
F_early_long_AFR$`maternal rank` <- scale(F_early_long_AFR$`maternal rank`)

AFR_cox_model <- coxph(AFR_object ~  `maternal rank` + `Ancylostoma egg load` + `f-IgA` + `f-GCM` , data = F_early_long_AFR)
summary(AFR_cox_model)


####### Hazard ratios plot - Fig.2
#select data
coefs <- AFR_cox_model$coefficients
var_covar <- vcov(AFR_cox_model)

# Calculate confidence intervals
plot_data <- data.frame(
  term = names(coefs),
  HR = exp(coefs),
  lower = exp(coefs - 1.96 * sqrt(diag(var_covar))),
  upper = exp(coefs + 1.96 * sqrt(diag(var_covar)))
)

# force order of predictors
model_order <- names(coef(AFR_cox_model))
plot_data$term <- factor(plot_data$term, levels = model_order)


#plot
forest_plot <- ggplot(plot_data, aes(x = HR, y = term)) +
  geom_point(size = 3, shape = 15) +  # Square points
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "orange") +
  scale_x_log10(limits = c(0.5, 4)) +
  labs(title = NULL, x = "Hazard Ratio", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size=14)
  )

forest_plot <- forest_plot +
  scale_y_discrete(labels = function(x) gsub("`", "", x))

forest_plot <- forest_plot +
  geom_text(
    aes(label = sprintf("%.2f (%.2f - %.2f)", HR, lower, upper), x = HR, y = term),
    nudge_y = -0.25,
    hjust = 0.5,
    vjust = 1,
    size = 4
  )
forest_plot

#save file
ggsave("HR_AFR.png",
       plot = forest_plot,
       width = 210,
       height = 148.5,
       units = "mm",
       dpi = 300)

####### LRT from anova
## fit full model and variants
m1 <- coxph(AFR_object ~  `maternal rank` + `Ancylostoma egg load` + `f-IgA` + `f-GCM` , data = F_early_long_AFR)
m0 <- update(m1, . ~ 1)
m2 <- coxph(AFR_object ~  `Ancylostoma egg load`  + `f-IgA` +  `f-GCM` , data = F_early_long_AFR)
m3 <- coxph(AFR_object ~  `maternal rank` +  `f-IgA` +  `f-GCM` , data = F_early_long_AFR)
m4 <- coxph(AFR_object ~  `maternal rank` + `Ancylostoma egg load` + `f-GCM` , data = F_early_long_AFR)
m5 <- coxph(AFR_object ~  `maternal rank` + `Ancylostoma egg load` +  `f-IgA`, data = F_early_long_AFR)


## test
anova(m0,m1) # intercept
anova(m2,m1) #maternal rank
anova(m3,m1) #Ancylostoma egg load
anova(m4,m1) #f-IgA
anova(m5,m1) # f-GCM

#####################################################
## Confirmation 6m

# restrict dataset to observations within 0.5 years of living
F_early_long_AFR <- subset(F_early_long, Survival_Ad == "Yes" )
F_early_long_AFR_6m <- subset(F_early_long_AFR,F_early_long_AFR$sample_days<= 183)
dim(F_early_long_AFR_6m) #34 obs
length(unique(F_early_long_AFR_6m$ID)) # 34 female individuals
sum(duplicated(F_early_long_AFR_6m)) # 0  duplicates

#create event variable:  reproduced or not
F_early_long_AFR_6m <- F_early_long_AFR_6m %>%
  mutate(Reproduced = case_when(AFR >= 0  ~ 1,
                                is.na(AFR) ~ 0))

#create variables for time to event variable (YEAR)
#females with AFR -> AFR
#females never reproduced and dead -> age of death
#female never reproduced still alive -> current date of dataset - birthdate
current_date = as.Date("2023-03-23", format = "%Y-%m-%d")
library(lubridate)
int <- interval(ymd(F_early_long_AFR_6m$birthdate), ymd(F_early_long_AFR_6m$deathdate))
F_early_long_AFR_6m$age_death <- time_length(int, "year")

int1 <- interval(ymd(F_early_long_AFR_6m$birthdate), ymd(current_date))
F_early_long_AFR_6m$current_age <- time_length(int1, "year")


# create time to event variable
F_early_long_AFR_6m<- F_early_long_AFR_6m%>%
  mutate(tafr = case_when(AFR >=0 ~ AFR,
                          is.na(AFR) & age_death >= 0 ~ age_death,
                          is.na(AFR) & is.na(age_death) ~ current_age))


# Fit a Cox proportional hazards model
# scale variables #
F_early_long_AFR_6m$`Ancylostoma egg load` <- scale(F_early_long_AFR_6m$`Ancylostoma egg load`)
F_early_long_AFR_6m$`f-IgA` <- scale(F_early_long_AFR_6m$`f-IgA`)
F_early_long_AFR_6m$`f-GCM` <- scale(F_early_long_AFR_6m$`f-GCM`)
F_early_long_AFR_6m$`maternal rank 6m` <- scale(F_early_long_AFR_6m$`maternal rank 6m`)


AFR_object_6m <- Surv(time = F_early_long_AFR_6m$tafr, event = F_early_long_AFR_6m$Reproduced)
AFR_object_6m
AFR_model_6m <- coxph(AFR_object_6m ~ `maternal rank 6m` +  `Ancylostoma egg load` +  `f-IgA` + `f-GCM` , data = F_early_long_AFR_6m)
summary(AFR_model_6m)



####### LRT from anova
## fit full model and variants
m1 <- coxph(AFR_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` + `f-IgA` + `f-GCM` , data = F_early_long_AFR_6m)
m0 <- update(m1, . ~ 1)
m2 <- coxph(AFR_object_6m ~  `Ancylostoma egg load` + `f-IgA` + `f-GCM`, data =F_early_long_AFR_6m)
m3 <- coxph(AFR_object_6m ~  `maternal rank 6m` +  `f-IgA` +  `f-GCM`, data = F_early_long_AFR_6m)
m4 <- coxph(AFR_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` +  `f-GCM` , data = F_early_long_AFR_6m)
m5 <- coxph(AFR_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` +  `f-IgA` , data = F_early_long_AFR_6m)


## test
anova(m0,m1) # intercept
anova(m2,m1) #maternal rank
anova(m3,m1) #Ancylostoma egg load
anova(m4,m1) #f-IgA
anova(m5,m1) #f-GCM



##################################################################################
### 1.3 Longevity                                               ####
##################################################################################
#removal of Ancylostoma egg load outlier
F_early_longv <- subset(F_early_long,`Ancylostoma egg load` < 21000 )

#create event variable: died or not
F_early_longv$Longevity[is.na(F_early_longv$Longevity)] <- "-"
F_early_longv<- F_early_longv %>%
  mutate(Died = case_when(Longevity >= 0 ~ 1,
                          Longevity == "-" ~ 0))
F_early_longv$Died <- as.numeric(F_early_longv$Died)


#create variable time to event (YEAR)
F_early_longv$Longevity <- as.numeric(F_early_longv$Longevity)
F_early_longv$sample_years <- as.numeric(F_early_longv$sample_years)
#females with death date -> age_death
#females without death date  -> current_age
current_date = as.Date("2023-03-23", format = "%Y-%m-%d")
int <- interval(ymd(F_early_longv$birthdate), ymd(F_early_longv$deathdate))
F_early_longv$age_death_years <- time_length(int, "year")

int1 <- interval(ymd(F_early_longv$birthdate), ymd(current_date))
F_early_longv$current_age <- time_length(int1, "year")

F_early_longv$age_death_years <- as.numeric(F_early_longv$age_death_years)
F_early_longv$current_age <- as.numeric(F_early_longv$current_age)

# create time to event variable
F_early_longv<- F_early_longv %>%
  mutate(ttdeath = case_when(deathdate >=0 ~ age_death_years,
                             TRUE ~ current_age))

F_early_longv$ttdeath <- as.numeric(F_early_longv$ttdeath)


# Fit the Cox - surv_object
Long_object <- Surv(time = F_early_longv$ttdeath, event = F_early_longv$Died)
Long_object

# Fit a Cox proportional hazards model with normalized data
#scale to mean 0 and sd 1
F_early_longv$`Ancylostoma egg load`<- scale(F_early_longv$`Ancylostoma egg load`)
F_early_longv$`f-mucin`<- scale(F_early_longv$`f-mucin`)
F_early_longv$`f-IgA`<- scale(F_early_longv$`f-IgA`)
F_early_longv$Polyparasitism <-scale(F_early_longv$Polyparasitism)
F_early_longv$`f-GCM` <- scale(F_early_longv$`f-GCM`)
F_early_longv$`maternal rank` <- scale(F_early_longv$`maternal rank`)

cox_longevity <- coxph(Long_object ~ `maternal rank` + `Ancylostoma egg load` + Polyparasitism + `f-IgA` + `f-mucin` + `f-GCM`, data = F_early_longv)
summary(cox_longevity)


####### Hazard ratios plot - Fig. 3
# select data
coefs <- cox_longevity$coefficients
var_covar <- vcov(cox_longevity)

# Calculate confidence intervals
plot_data <- data.frame(
  term = names(coefs),
  HR = exp(coefs),
  lower = exp(coefs - 1.96 * sqrt(diag(var_covar))),  # Lower CI
  upper = exp(coefs + 1.96 * sqrt(diag(var_covar)))   # Upper CI
)

#order predictors
model_order <- names(coef(cox_longevity))
plot_data$term <- factor(plot_data$term, levels = model_order)


#plot Hazard ratios
forest_plot <- ggplot(plot_data, aes(x = HR, y = term)) +
  geom_point(aes(color = ifelse(term == "`maternal rank`", "darkorchid", "black")), size = 3, shape = 15) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = ifelse(term == "`maternal rank`", "darkorchid", "black")), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "orange") +
  scale_x_log10(limits = c(0.5, 4)) +
  labs(title = NULL, x = "Hazard Ratio", y = NULL) +
  theme_minimal() +
  scale_color_manual(values = c("darkorchid" = "darkorchid", "black" = "black")) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size=14),
    legend.position = "none"
  )

forest_plot <- forest_plot +
  scale_y_discrete(labels = function(x) gsub("`", "", x))

forest_plot <- forest_plot +
  geom_text(
    aes(label = sprintf("%.2f (%.2f - %.2f)", HR, lower, upper), x = HR, y = term),
    nudge_y = -0.25,
    hjust = 0.5,
    vjust = 1,
    size = 4
  )
forest_plot

#save file
ggsave("HR_LV.png",
       plot = forest_plot,
       width = 210,
       height = 148.5,
       units = "mm",
       dpi = 300)

####### LRT from anova
## fit full model and variants
m1 <- coxph(Long_object ~  `maternal rank` + `Ancylostoma egg load` +  Polyparasitism +  `f-IgA` + `f-mucin` + `f-GCM` , data = F_early_longv)
m0 <- update(m1, . ~ 1)
m2 <- coxph(Long_object ~  `Ancylostoma egg load` +  Polyparasitism +  `f-IgA` + `f-mucin` + `f-GCM` , data = F_early_longv)
m3 <- coxph(Long_object ~  `maternal rank` +  Polyparasitism +  `f-IgA` + `f-mucin` + `f-GCM` , data = F_early_longv)
m4 <- coxph(Long_object ~  `maternal rank` + `Ancylostoma egg load` +  `f-IgA` + `f-mucin` + `f-GCM` , data = F_early_longv)
m5 <- coxph(Long_object ~  `maternal rank` + `Ancylostoma egg load` + Polyparasitism +  `f-mucin` + `f-GCM` , data = F_early_longv)
m6 <- coxph(Long_object ~  `maternal rank` + `Ancylostoma egg load` +  Polyparasitism +  `f-IgA` + `f-GCM` , data = F_early_longv)
m7 <- coxph(Long_object ~  `maternal rank` + `Ancylostoma egg load` +  Polyparasitism +  `f-IgA` + `f-mucin`, data = F_early_longv)


## test
anova(m0,m1) # intercept
anova(m2,m1) #maternal rank
anova(m3,m1) #Ancylostoma egg load
anova(m4,m1) #Polyparasitsm
anova(m5,m1) #f-IgA
anova(m6,m1) #f-mucin
anova(m7,m1) #f-GCM

#############plot significant effect on adjusted curves -  Fig. 4
F_early_longv <- subset(F_early_long,`Ancylostoma egg load` < 21000 )
F_early_longv_plot <- F_early_longv %>%
  rename (maternal_rank = `maternal rank`,
          Ancylostoma_egg_load = `Ancylostoma egg load`,
          f_IgA = `f-IgA`, f_mucin = `f-mucin`, f_GCM = `f-GCM`)

F_early_longv_plot$maternal_rank_group <- factor(ifelse(F_early_longv_plot$maternal_rank < 0,
                                                        "Low Rank", "High Rank"), levels = c("Low Rank", "High Rank"))

#create event variable: died or not
F_early_longv_plot$Longevity[is.na(F_early_longv_plot$Longevity)] <- "-"
F_early_longv_plot<- F_early_longv_plot %>%
  mutate(Died = case_when(Longevity >= 0 ~ 1,
                          Longevity == "-" ~ 0))
F_early_longv_plot$Died <- as.numeric(F_early_longv_plot$Died)


#create variable time to event (YEAR)
F_early_longv_plot$Longevity <- as.numeric(F_early_longv_plot$Longevity)
F_early_longv_plot$sample_years <- as.numeric(F_early_longv_plot$sample_years)
#females with death date -> age_death - sample_age(years)
#females without death date  -> current_age - sample_age(years)
#removed sample age because we are taking one sample per ID
current_date = as.Date("2023-03-23", format = "%Y-%m-%d")
int <- interval(ymd(F_early_longv_plot$birthdate), ymd(F_early_longv_plot$deathdate))
F_early_longv_plot$age_death_years <- time_length(int, "year")

int1 <- interval(ymd(F_early_longv_plot$birthdate), ymd(current_date))
F_early_longv_plot$current_age <- time_length(int1, "year")

F_early_longv_plot$age_death_years <- as.numeric(F_early_longv_plot$age_death_years)
F_early_longv_plot$current_age <- as.numeric(F_early_longv_plot$current_age)

# create time to event variable
F_early_longv_plot<- F_early_longv_plot %>%
  mutate(ttdeath = case_when(deathdate >=0 ~ age_death_years,
                             TRUE ~ current_age))

F_early_longv_plot$ttdeath <- as.numeric(F_early_longv_plot$ttdeath)


# Fit the Cox - surv_object
Long_object <- Surv(time = F_early_longv_plot$ttdeath, event = F_early_longv_plot$Died)
Long_object

cox_longevity <- coxph(Long_object ~ maternal_rank_group + Ancylostoma_egg_load + Polyparasitism + f_IgA + f_mucin + f_GCM, data = F_early_longv_plot)
summary(cox_longevity)

# Create a simplified data frame for plotting
F_early_longv_plot2 <- data.frame(maternal_rank_group = c("Low Rank", "High Rank"))

# Add other covariates from the model, setting them to their meadian value
F_early_longv_plot2$f_IgA <- median(F_early_longv_plot$f_IgA)
F_early_longv_plot2$f_mucin <- median(F_early_longv_plot$f_mucin)
F_early_longv_plot2$f_GCM <- median(F_early_longv_plot$f_GCM)
F_early_longv_plot2$Ancylostoma_egg_load <- median(F_early_longv_plot$Ancylostoma_egg_load)
F_early_longv_plot2$Polyparasitism <- median(F_early_longv_plot$Polyparasitism)

#use surfit to create curves
surv_preds <- survfit(cox_longevity, newdata = F_early_longv_plot2)
surv_data_plot <- tidy(surv_preds)
names(surv_data_plot)

surv_data_long <- surv_data_plot %>%
  pivot_longer(
    cols = starts_with("estimate"),
    names_to = "group",
    values_to = "estimate"
  ) %>%
  mutate(
    group = sub("estimate.", "", group)
  )

surv_data_long <- surv_data_long %>%
  left_join(
    surv_data_plot %>%
      pivot_longer(
        cols = starts_with("conf.high"),
        names_to = "group_h",
        values_to = "conf.high"
      ) %>%
      select(-group_h),
    by = "time"
  ) %>%
  left_join(
    surv_data_plot %>%
      pivot_longer(
        cols = starts_with("conf.low"),
        names_to = "group_l",
        values_to = "conf.low"
      ) %>%
      select(-group_l),
    by = "time"
  ) %>%
  mutate(
    group = as.factor(group)
  )

#plot
Adj_curve <- ggplot(surv_data_long, aes(x = time, y = estimate, color = group)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (Years)",
    y = "Adjusted Survival probability",
    color = "maternal rank"
  ) +
  scale_y_continuous(labels = percent) +
  scale_color_manual(
    values = c("1" = "#D8BFD8", "2" = "darkorchid"),
    labels = c("1" = "Low Ranking", "2" = "High Ranking")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

#save file
ggsave("Adj_maternal.png",
       plot = Adj_curve,
       width = 210,
       height = 148.5,
       units = "mm",
       dpi = 300)


#####################################
## confirmation 6m ##
# restrict dataset to observations within 0.5 years of living
F_early_longv_6m <- subset(F_early_long,F_early_long$sample_days<= 183)
dim(F_early_longv_6m) #51obs
length(unique(F_early_longv_6m$ID)) # 51 female individuals
sum(duplicated(F_early_longv_6m)) # 0  duplicates

F_early_longv_6m$`Ancylostoma egg load`<- scale(F_early_longv_6m$`Ancylostoma egg load`)
F_early_longv_6m$`f-mucin`<- scale(F_early_longv_6m$`f-mucin`)
F_early_longv_6m$`f-IgA`<- scale(F_early_longv_6m$`f-IgA`)
F_early_longv_6m$Polyparasitism <-scale(F_early_longv_6m$Polyparasitism)
F_early_longv_6m$`f-GCM`<- scale(F_early_longv_6m$`f-GCM`)
F_early_longv_6m$`maternal rank 6m` <- scale(F_early_longv_6m$`maternal rank 6m`)

#create event variable: died or not
F_early_longv_6m$Longevity[is.na(F_early_longv_6m$Longevity)] <- "-"
F_early_longv_6m<- F_early_longv_6m %>%
  mutate(Died = case_when(Longevity >= 0 ~ 1,
                          Longevity == "-" ~ 0))
F_early_longv_6m$Died <- as.numeric(F_early_longv_6m$Died)


#create variable time to event (YEAR)
F_early_longv_6m$Longevity <- as.numeric(F_early_longv_6m$Longevity)
F_early_longv_6m$sample_years <- as.numeric(F_early_longv_6m$sample_years)
#females with death date -> age_death
#females without death date  -> current_age

current_date = as.Date("2023-03-23", format = "%Y-%m-%d")
int <- interval(ymd(F_early_longv_6m$birthdate), ymd(F_early_longv_6m$deathdate))
F_early_longv_6m$age_death_years <- time_length(int, "year")

int1 <- interval(ymd(F_early_longv_6m$birthdate), ymd(current_date))
F_early_longv_6m$current_age <- time_length(int1, "year")

F_early_longv_6m$age_death_years <- as.numeric(F_early_longv_6m$age_death_years)
F_early_longv_6m$current_age <- as.numeric(F_early_longv_6m$current_age)

# create time to event variable
F_early_longv_6m<- F_early_longv_6m %>%
  mutate(ttdeath = case_when(deathdate >=0 ~ age_death_years,
                             TRUE ~ current_age))

F_early_longv_6m$ttdeath <- as.numeric(F_early_longv_6m$ttdeath)

# Fit the Cox - surv_object
Long_object_6m <- Surv(time = F_early_longv_6m$ttdeath, event = F_early_longv_6m$Died)
Long_object_6m

# Fit a Cox proportional hazards model

cox_longevity_6m <- coxph(Long_object_6m ~ `maternal rank 6m` + `Ancylostoma egg load` + Polyparasitism + `f-IgA` + `f-mucin` + `f-GCM`,data = F_early_longv_6m)
summary(cox_longevity_6m)

## LRT from anova
## fit full model and variants
m1 <- coxph(Long_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` +  Polyparasitism +  `f-IgA` + `f-mucin` + `f-GCM` , data = F_early_longv_6m)
m0 <- update(m1, . ~ 1)
m2 <- coxph(Long_object_6m ~  `Ancylostoma egg load` +  Polyparasitism +  `f-IgA` + `f-mucin` + `f-GCM` , data = F_early_longv_6m)
m3 <- coxph(Long_object_6m ~  `maternal rank 6m` +  Polyparasitism +  `f-IgA` + `f-mucin` + `f-GCM` , data = F_early_longv_6m)
m4 <- coxph(Long_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` +  `f-IgA` + `f-mucin` + `f-GCM`, data = F_early_longv_6m)
m5 <- coxph(Long_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` + Polyparasitism +  `f-mucin` + `f-GCM`, data = F_early_longv_6m)
m6 <- coxph(Long_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` +  Polyparasitism +  `f-IgA` + `f-GCM`, data = F_early_longv_6m)
m7 <- coxph(Long_object_6m ~  `maternal rank 6m` + `Ancylostoma egg load` +  Polyparasitism +  `f-IgA` + `f-mucin`, data = F_early_longv_6m)


## test
anova(m0,m1) #intercept
anova(m2,m1) #maternal rank
anova(m3,m1) #Ancylostoma egg load
anova(m4,m1) #Polyparasitism
anova(m5,m1) #f-IgA
anova(m6,m1) #f-mucin
anova(m7,m1) #f-GCM



sessionInfo()
RStudio.Version()
citation()
