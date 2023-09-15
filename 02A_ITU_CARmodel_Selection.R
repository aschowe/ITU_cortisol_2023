# ITU: CAR Model Selection and Baseline Description ############
# 18.11.22
# alicia_schowe@pysch.mpg.de 

# Data --------------------------------------------------------------------
load("Data/CAR.df.Rdata")

# libraries ---------------------------------------------------------------
library(nlme)
library(data.table)
library(APAstyler) #remotes::install_github("darrellpenta/APAstyler")


# Functions ---------------------------------------------------------------
source("ITU_lmm_helper_functions.R")


# Baseline Model: Random effects and variance-covariance structure ---------
#Including timing factors 
CAR_Form <- lg_cort_concentration_ug_l ~ 
  hours_since_waking +
  I(hours_since_waking^2) + 
  pregstage + 
  waking_time_cent + 
  season +
  hours_since_waking:pregstage +
  hours_since_waking:season + 
  hours_since_waking:waking_time_cent + 
  pregstage:season + 
  pregstage:waking_time_cent + 
  hours_since_waking:pregstage:season + 
  hours_since_waking:pregstage:waking_time_cent
  
## 1. 3- vs 4-level modeling? 
## Random Intercept only model 

### 3 levels 
M.CAR.RI.3level <- lme(CAR_Form,
                       random = list(~1|ELISA_analysis_plate,
                                     ~ 1|participantID),
                       method = "REML",
                       data = CAR.df,
                       na.action = "na.exclude",
                       control= lmeControl(msMaxIter = 200))

plot(M.CAR.RI.3level)
qqnorm(M.CAR.RI.3level, abline = c(0,1))

### 4levels 
M.CAR.RI.4level <- update(M.CAR.RI.3level,
                       random = list(~1|ELISA_analysis_plate,
                                     ~ 1|participantID,
                                     ~ 1|pregstage))
plot(M.CAR.RI.4level)
qqnorm(M.CAR.RI.4level, abline = c(0,1))

CAR_RES <- anova(M.CAR.RI.3level, M.CAR.RI.4level) #4 level formula has sig. better fit 

## 2. Does adding a random slope for CAR improve model fit? RI + RS model 
#i.e., allowing morning slope/usually=CAR to vary across participants 
M.CAR.RS <- update(M.CAR.RI.3level, random = list(~1|ELISA_analysis_plate,
                                                  ~hours_since_waking|participantID,
                                                  ~1|pregstage))
res.RS <- anova(M.CAR.RI.4level, M.CAR.RS) #Yes, it does improve fit 
CAR_RES <- rbind.data.frame(CAR_RES, res.RS)


# 3. Does adding heterogenous WS variance for the effect of time improve fit?
# i.e., controlling for differential data variation at different times of the day (morning sample have greater spread)
M.CAR.hetero <- update(M.CAR.RS,
                    random = list(~1|ELISA_analysis_plate,
                                  ~hours_since_waking|participantID,
                                  ~1|pregstage),
                    weights = varExp(form =  ~hours_since_waking))

plot(M.CAR.hetero)
res.hetero <- anova(M.CAR.RS, M.CAR.hetero) #Does not sig. improve fit!
CAR_RES <- rbind.data.frame(CAR_RES, res.hetero)

# 4. Does adding auto-correlation improve fit? 
# i.e., controlling for stronger correlation when samples are closer together in time 
M.CAR.auto <- update(M.CAR.RS,
                  random = list(~1|ELISA_analysis_plate,
                                ~hours_since_waking|participantID,
                                ~1|pregstage),
                  correlation = corCAR1(form =  ~ (hours_since_waking)))
plot(M.CAR.auto)

res.auto <- as.data.frame(anova(M.CAR.RS, M.CAR.auto)) #Does improve fit!
CAR_RES_Mselect <- rbind.data.frame(anova(M.CAR.RI.3level, 
                         M.CAR.RI.4level, 
                         M.CAR.RS,
                         M.CAR.hetero),
                         res.auto)[-5,]  #remove doubled row

CAR_RES_Mselect[5,"Model"] <- c("5")
CAR_RES_Mselect$Test <- as.character(CAR_RES_Mselect$Test)
CAR_RES_Mselect[5,"Test"] <- c("3 vs 5")

## Model selection results ---------------------------
save(CAR_RES_Mselect, file="results/CAR_RES_Mselect.Rdata")

# Basic Model: within-subject time fixed effects -------------------------------
M.CAR_baseline <- update(M.CAR.auto, method = "ML")
chi2_CAR <- as.data.frame(car::Anova(M.CAR_baseline))
chi2_CAR$Chisq <- round(chi2_CAR$Chisq, digits = 1)
chi2_CAR$`Pr(>Chisq)` <- round(chi2_CAR$`Pr(>Chisq)`, digits = 3)
chi2_CAR$`Pr(>Chisq)`[chi2_CAR$`Pr(>Chisq)` < 0.001] <- rep("<.001", length(chi2_CAR$`Pr(>Chisq)`[chi2_CAR$`Pr(>Chisq)` < 0.001]))

save(chi2_CAR, file = "results/Tables/CAR_chi2_Fullbasic.Rdata")
writexl::write_xlsx(chi2_CAR, "results/Tables/CAR_chi2_Fullbasic.xlsx")

# take out n.s. 3-way interactions 
M.CAR_baseline_3way <- update(M.CAR_baseline,
                               .~. - hours_since_waking:pregstage:season -
                                hours_since_waking:pregstage:waking_time_cent)
car::Anova(M.CAR_baseline_3way)

# take out n.s. 2 interaction effects, ie., pregstage:season, 
# pregstage:waking_time_cent, hours_since_waking:season, hours_since_waking:waking_time_cent
M.CAR_baseline_final <- update(M.CAR_baseline_3way,
                               .~. - pregstage:season -
                                 pregstage:waking_time_cent -
                                 hours_since_waking:season -
                                 hours_since_waking:waking_time_cent)

car::Anova(M.CAR_baseline_final)
summary(M.CAR_baseline_final)

# format output -------------------------------------------------------------
RES.CAR_Baseline <- model_res.wRandom(summary(M.CAR_baseline_final))
RES.CAR_Baseline <- mutate_all(RES.CAR_Baseline,funs(replace(., is.na(.), "")))
#RES.CAR_Baseline$Interpretation <- NULL
RES.CAR_Baseline$predictor[1:9] <- c("Intercept",
                                "Time since awakening",
                                "Time since awakening ^2",
                                "T2 (ref = T1)",
                                "T3 (ref = T1)",
                                "Time at awakening", 
                                "Season [winter]",
                                "Time since awakening x T2",
                                "Time since awakening x T3")

# ICC in Basic Model ----------------------------------------------------------
## add results to baseline model 

ICC_within_pregstage <- ICC_within(M.CAR_baseline_final)
ICC_across_preg <- ICC_across(M.CAR_baseline_final) 
ICC_df <- data.frame(matrix(ncol = 6, nrow = 3))
ICC_df[1,] <- c("ICC estimates", rep("", 5))
ICC_df[2,] <- c("ICC across occasions", ICC_across_preg, rep("", 3), "Estimated stability of cortisol levels of the same individual across mornings")
ICC_df[3,] <- c("ICC within occasions", ICC_within_pregstage, rep("", 3), "Estimated correlation of cortisol levels of the same individual within mornings")
names(ICC_df) <- names(RES.CAR_Baseline)
RES.CAR_Baseline <- rbind(RES.CAR_Baseline, ICC_df)
save(RES.CAR_Baseline, file="results/CAR_BasicM_RES.Rdata")
writexl::write_xlsx(RES.CAR_Baseline, "results/Tables/CAR_BasicM_RES_02052023.xlsx")

# relevel to compare to mid-pregnancy for in-text results
CAR.df$pregstage <- relevel(factor(CAR.df$pregstage, 
                           ordered=F), 
                           ref = "T2")
M.CAR_baselineT2 <- update(M.CAR_baseline_final, 
                          data = CAR.df)
RES.CAR_BaselineII <- summary(M.CAR_baselineT2)

pregstageIIvsIII <- RES.CAR_BaselineII$tTable[5,]
slopeIIvsIII <- RES.CAR_BaselineII$tTable[9,]

