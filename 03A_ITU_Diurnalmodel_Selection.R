# ITU: Diurnal Model Selection and Baseline results ######################
# 27.11.22
# alicia_schowe@psych.mpg.de 

# Data  -----------------------------------------------------------------------
load("Data/diurnal.df.Rdata")

# libraries  -------------------------------------------------------------------
library(nlme)
library(data.table)
library(APAstyler)

# Functions  -------------------------------------------------------------------
source("ITU_lmm_helper_functions.R")

# Baseline Model: Random effects and variance-covariance structure -------------
DIURNAL_Form <- lg_cort_concentration_ug_l ~ 
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
  pregstage:season:hours_since_waking + 
  pregstage:waking_time_cent:hours_since_waking

## 1. Random Intercept only model 
M.DIUR.RI1 <- lme(DIURNAL_Form,
                  random = list(~1|ELISA_analysis_plate,
                                ~ 1|participantID),
                  method = "REML",
                  data = diurnal.df,
                  na.action = "na.omit",
                  control= lmeControl(msMaxIter = 200))
#plot(M.DIUR.RI1)
#qqnorm(M.DIUR.RI1, abline = c(0,1))

## Does adding a third level improve model fit? 
M.DIUR.RI2 <- update(M.DIUR.RI1, random = list(~1|ELISA_analysis_plate,
                                               ~1|participantID,
                                               ~1|pregstage))

RES <- anova(M.DIUR.RI1,M.DIUR.RI2) #Yes, improved model fit 

## 2. Does adding a random slope for time since waking improve model fit? 
#i.e., allowing the diurnal decline to vary across participants 
# M.DIUR.RS1 <- update(M.DIUR.RI2, random = list(~1|ELISA_analysis_plate,
#                                                 ~hours_since_waking|participantID,
#                                                 ~1|pregstage)) #singular convergence


# 3. Does adding heterogenous WS variance for the effect of time improve fit?
M.DIUR.hetero <- lme(DIURNAL_Form,
                     random = list(~1|ELISA_analysis_plate,
                                   ~1|participantID,
                                   ~1|pregstage),
                     weights = varExp(form =  ~hours_since_waking),
                     method = "REML",
                     data = diurnal.df,
                     na.action = "na.omit",
                     control= lmeControl(msMaxIter = 200))
anova(M.DIUR.RI2, M.DIUR.hetero) #Yes, it does improve fit!
#plot(M.DIUR.hetero)

### 4. Can the slope for cortisol decline now be added without convergence issues 
#Does adding heterogenous WS variance for the effect of time improve fit?
M.DIUR.hetero2 <- lme(DIURNAL_Form,
                      random = list(~1|ELISA_analysis_plate,
                                    ~(hours_since_waking)|participantID,
                                    ~1|pregstage),
                      weights = varExp(form =  ~hours_since_waking),
                      method = "REML",
                      data = diurnal.df,
                      na.action = "na.omit",
                      control= lmeControl(msMaxIter = 200))
anova(M.DIUR.hetero, M.DIUR.hetero2) #Yes, it does improve fit!
#plot(M.DIUR.hetero)

# 5. Does adding autocorrelation for WS variance for the effect of time improve fit?
M.DIUR.hetero_auto <- lme(DIURNAL_Form,
                          random = list(~1|ELISA_analysis_plate,
                                        ~(hours_since_waking)|participantID,
                                        ~1|pregstage),
                          weights = varExp(form =  ~hours_since_waking|pregstage),
                          correlation = corCAR1(form =  ~ (hours_since_waking)),
                          method = "REML",
                          data = diurnal.df,
                          na.action = "na.omit",
                          control= lmeControl(msMaxIter = 200))
plot(M.DIUR.hetero_auto)
qqnorm(M.DIUR.hetero_auto, abline = c(0,1))

## Baseline model selection results -------------------------------------------
RES_DIUR_Mselect <- anova(M.DIUR.RI1, M.DIUR.RI2,M.DIUR.hetero, M.DIUR.hetero2,M.DIUR.hetero_auto) #Yes, it does improve fit!
save(RES_DIUR_Mselect, file="results/DIUR_RES_Mselect.Rdata")

# Baseline Model: Fixed Component -------------------------------------------
M.DIUR_baseline <- update(M.DIUR.hetero_auto, method = "ML")
chi2_DIUR <- as.data.frame(car::Anova(M.DIUR_baseline))
chi2_DIUR$Chisq <- round(chi2_DIUR$Chisq, digits = 1)
chi2_DIUR$`Pr(>Chisq)` <- round(chi2_DIUR$`Pr(>Chisq)`, digits = 3)
chi2_DIUR$`Pr(>Chisq)`[chi2_DIUR$`Pr(>Chisq)` < 0.001] <- rep("<.001", length(chi2_DIUR$`Pr(>Chisq)`[chi2_DIUR$`Pr(>Chisq)` < 0.001]))

save(chi2_DIUR, file = "results/DIUR_chi2_Fullbasic.Rdata")
writexl::write_xlsx(chi2_DIUR, "results/Tables/chi2_DIUR_Fullbasic.xlsx")


#### refit after excluding n.s. 3-way interactions 
M.DIUR_baseline_3way <- update(M.DIUR_baseline, .~. -(hours_since_waking:pregstage:season))
car::Anova(M.DIUR_baseline_3way, type="III")

#### refit after excluding n.s. 2-way interactions
M.DIUR_baseline_final <- update(M.DIUR_baseline_3way, .~. -(hours_since_waking:season))
car::Anova(M.DIUR_baseline_final, type="III")


# format output -------------------------------------------------------------
car::Anova(M.DIUR_baseline_final)
RES_DIUR_Baseline <- model_res.wRandom(summary(M.DIUR_baseline_final))
RES_DIUR_Baseline <- mutate_all(RES_DIUR_Baseline,funs(replace(., is.na(.), "")))

#### Null
ICC_day <- ICC_within(M.DIUR_baseline_final)
ICC_preg <- ICC_across(M.DIUR_baseline_final)

## add results to baseline model 
ICC_df <- data.frame(matrix(ncol = 6, nrow = 3))
ICC_df[1,] <- c("ICC estimates", rep("", 5))
ICC_df[2,] <- c("ICC across occasions", ICC_preg, rep("", 3), "Estimated stability of cortisol levels of the same individual across pregnancy")
ICC_df[3,] <- c("ICC within occasions", ICC_day, rep("", 3), "Estimated correlation of cortisol levels of the same individual within sampling days")
names(ICC_df) <- names(RES_DIUR_Baseline)
RES_DIUR_Baseline <- rbind(RES_DIUR_Baseline, ICC_df)
#save(RES_DIUR_Baseline, file="results/DIUR_BasicM_RES.Rdata")
writexl::write_xlsx(RES_DIUR_Baseline, "results/Tables/RES_DIUR_Baselinebasic.xlsx")

# Relevel baseline model to compare effects to mid-pregnancy -----------------
diurnal.df$pregstage <- factor(diurnal.df$pregstage, 
                               ordered = F)
M.Baseline_DIUR_cort2_II <- update(M.DIUR_baseline_final, 
                              lg_cort_concentration_ug_l ~ 
                                hours_since_waking +
                                I(hours_since_waking^2) + 
                                relevel(pregstage, ref = "T2") + 
                                waking_time_cent + 
                                season +
                                hours_since_waking:relevel(pregstage, ref = "T2") +
                                hours_since_waking:waking_time_cent +
                                relevel(pregstage, ref = "T2"):waking_time_cent + 
                                relevel(pregstage, ref = "T2"):waking_time_cent:hours_since_waking)
DIUR_baselineII <- summary(M.Baseline_DIUR_cort2_II)

pregstageIIvsIII <- DIUR_baselineII$tTable[5,]
slopeIIvsIII <- DIUR_baselineII$tTable[9,]
waking_timeIIvsIII <- DIUR_baselineII$tTable[14,]



