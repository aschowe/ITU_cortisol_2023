# 04 Stratified analyses for BMI 

# setup ==================================================================
library(nlme)
library(data.table)
library(APAstyler) #remotes::install_github("darrellpenta/APAstyler")
library(naniar)
library(dplyr)
library(ggplot2)

load("Data/diurnal.df.Rdata")
load("Data/CAR.df.Rdata")

covariates <- c("caseVScontrol",
                "Parity",
                "Maternal_Age_Years",
                "Maternal_Education",
                "Maternal_Smoking_During_Pregnancy",
                "Maternal_Mental_Disorders")

source("ITU_lmm_helper_functions.R")

# =============================================================================
  
# CAR model ===========================================================
M1_trait_formula <- 
  as.formula(lg_cort_concentration_ug_l ~ 
               hours_since_waking + 
               I(hours_since_waking^2) + 
               pregstage + 
               waking_time_cent + 
               season + 
               hours_since_waking:pregstage) 

M.CAR.full <- lme(M1_trait_formula,
                  random = list(~1|ELISA_analysis_plate,
                                ~hours_since_waking|participantID,
                                ~1|pregstage),
                  correlation = corCAR1(form =  ~ (hours_since_waking)),
                  na.action = "na.exclude",
                  method = "ML",
                  data = CAR.df)


## Stratified Statistics: BMI < 25 ----------------------
CAR_BMI_low <- subset(CAR.df, Maternal_Body_Mass_Index_in_Early_Pregnancy < 25)
M1 <- update(M.CAR.full, data = CAR_BMI_low)
RES.CAR_final_M1 <- model_res(summary(M1))
RES.CAR_final_M1 <- mutate_all(RES.CAR_final_M1,funs(replace(., is.na(.), "")))
rownames(RES.CAR_final_M1)[2:3] <- c("Time", "Timeˆ2")
kable(RES.CAR_final_M1[2:3,c(1:4)],
      digits = 3,
      caption = "Time slope for BMI < 25")


M2_form <- paste0(".~.", "+", paste(covariates, collapse = "+")) 
M2 <- update(M.CAR.full, M2_form, data = CAR_BMI_low)
RES.CAR_final_M2 <- model_res(summary(M2))
RES.CAR_final_M2 <- mutate_all(RES.CAR_final_M2,funs(replace(., is.na(.), "")))
rownames(RES.CAR_final_M2)[2:3] <- c("Time", "Timeˆ2")
RES.CAR_final_M2

##  Stratified Statistics BMI >= 25 ------------------------
CAR_BMI_high <- subset(CAR.df, Maternal_Body_Mass_Index_in_Early_Pregnancy >= 25)
M1 <- update(M.CAR.full, data = CAR_BMI_high)
RES.CAR_final_M1 <- model_res(summary(M1))
RES.CAR_final_M1 <- mutate_all(RES.CAR_final_M1,funs(replace(., is.na(.), "")))
rownames(RES.CAR_final_M1)[2:3] <- c("Time", "Timeˆ2")
kable(RES.CAR_final_M1[2:3,c(1:4)],
      digits = 3,
      caption = "Time slope for BMI > 25")


M2_form <- paste0(".~.", "+", paste(covariates, collapse = "+")) 
M2 <- update(M.CAR.full, M2_form, data = CAR_BMI_high)
RES.CAR_final_M2 <- model_res(summary(M2))
RES.CAR_final_M2 <- mutate_all(RES.CAR_final_M2,funs(replace(., is.na(.), "")))
rownames(RES.CAR_final_M2)[2:3] <- c("Time", "Timeˆ2")

# DIUR ======================================================================
M1_trait_formula <- 
  as.formula(lg_cort_concentration_ug_l ~ 
               hours_since_waking + 
               I(hours_since_waking^2) + 
               pregstage + 
               waking_time_cent + 
               season + 
               season:pregstage + 
               hours_since_waking:pregstage + 
               hours_since_waking:waking_time_cent + 
               hours_since_waking:waking_time_cent:pregstage)

M.DIUR.full <- lme(M1_trait_formula,
                   random = list(~1|ELISA_analysis_plate,
                                 ~1|participantID,
                                 ~1|pregstage),
                   weights = varExp(form =  ~hours_since_waking|pregstage),
                   correlation = corCAR1(form =  ~ (hours_since_waking)),
                   method = "ML",
                   data = diurnal.df,
                   na.action = "na.exclude",
                   control= lmeControl(msMaxIter = 200))

## BMI <= 25 ----------------------------
DCS_BMI_low <- subset(diurnal.df, Maternal_BMI_Early_Preg_normalVSover_obese == "BMI < 25")
M1 <- update(M.DIUR.full, data = DCS_BMI_low)
RES.DCS_final_M1 <- model_res(summary(M1))
RES.DCS_final_M1 <- mutate_all(RES.DCS_final_M1,funs(replace(., is.na(.), "")))
rownames(RES.DCS_final_M1)[10:11] <- c("Time x T2", "Time x T3")

kable(RES.DCS_final_M1[10:11,c(1:4)],
      digits = 3,
      caption = "Time x Pregnancy Stage for BMI < 25")

M2_form <- paste0(".~.", "+", paste(covariates, collapse = "+")) 
M2 <- update(M.DIUR.full, M2_form, data = DCS_BMI_low)
RES.DCS_final_M2 <- model_res(summary(M2))
RES.DCS_final_M2 <- mutate_all(RES.DCS_final_M2,funs(replace(., is.na(.), "")))
rownames(RES.DCS_final_M2)[17:18] <- c("Time x T2", "Time x T3")

## BMI >=25 -----------------------------
DCS_BMI_high <- subset(diurnal.df, Maternal_BMI_Early_Preg_normalVSover_obese == "BMI >= 25")
M1 <- update(M.DIUR.full, data = DCS_BMI_high)
RES.DCS_final_M1 <- model_res(summary(M1))
RES.DCS_final_M1 <- mutate_all(RES.DCS_final_M1,funs(replace(., is.na(.), "")))
rownames(RES.DCS_final_M1)[10:11] <- c("Time x T2", "Time x T3")

kable(RES.DCS_final_M1[10:11,c(1:4)],
      digits = 3,
      caption = "Time x Pregnancy Stage for BMI >= 25")

M2_form <- paste0(".~.", "+", paste(covariates, collapse = "+")) 
M2 <- update(M.DIUR.full, M2_form, data = DCS_BMI_high)
RES.DCS_final_M2 <- model_res(summary(M2))
RES.DCS_final_M2 <- mutate_all(RES.DCS_final_M2,funs(replace(., is.na(.), "")))
rownames(RES.DCS_final_M2)[17:18] <- c("Time x T2", "Time x T3")
