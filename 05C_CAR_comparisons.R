# ITU: Follow-up analyses of CAR findings 
# 31.03.2022
# alicia_schowe@psych.mpg.de 

# libraries ------------------------------------------------------------------
library(dplyr)
library(gtsummary)

# Data -----------------------------------------------------------------------
load("Data/ITU_combined_cortisol_dates_times_wide_format.Rdata")
wide_cort$peak_concentration <- apply(wide_cort[, c("lg_cort_concentration_ug_l_S1",
                                                    "lg_cort_concentration_ug_l_S2",
                                                    "lg_cort_concentration_ug_l_S3")], 
                                      1, max)
wide_cort$peak_timeslot <- NA
for(i in 1:nrow(wide_cort)){
  peak_con <- wide_cort$peak_concentration[[i]]
  if(!is.na(peak_con)){
  if(peak_con == wide_cort$lg_cort_concentration_ug_l_S1[[i]]){
    wide_cort$peak_timeslot[[i]] <- "S1" 
  }
  if(peak_con == wide_cort$lg_cort_concentration_ug_l_S2[[i]]){
    wide_cort$peak_timeslot[[i]] <- "S2" 
  }
  if(peak_con == wide_cort$lg_cort_concentration_ug_l_S3[[i]]){
    wide_cort$peak_timeslot[[i]] <- "S3" 
  }
  }
}

load("Data/register_dat.final.Rdata")
cort_phenos <- left_join(wide_cort, register_dat.final)

# (1) Compare women with and without CAR 
## Early Preg ----------------------------------------------------------------

demo <- as.data.frame(cort_phenos[cort_phenos$pregstage == "I" & cort_phenos$CAR_trend %in% c("positive", 
                                                                                              "negative"),]) %>% 
  select(c("CAR_trend", 
           "CAR_iAUC", 
           "gestage_weeks",
           "waking_time", 
           "lg_waking_cort_concentration_ug_l",
           "peak_concentration",
           "peak_timeslot",
           "season",
           "caseVScontrol",
           "Maternal_Age_Years",
           "Parity",
           "Maternal_Education",
           "Maternal_Smoking_During_Pregnancy",
           "Maternal_Mental_Disorders",
           "Maternal_Diabetes_Disorders_anyVSnone",
           "Maternal_Hypertensive_Disorders_anyVSnone",
           "Maternal_Body_Mass_Index_in_Early_Pregnancy",
           "weight_gain_during_preg",
           "relative_to_recomd_weight_gain"))

CAR.Tbl1 <- 
  tbl_summary(data = demo, 
              by = CAR_trend,
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"), 
              digits = list(all_continuous() ~ c(2, 2, 2)),           
              label = list(gestage_weeks ~ "Gestational Week at Cortisol Asssessment",
                           waking_time ~ "Time at Awakening",
                           lg_waking_cort_concentration_ug_l ~"Cortisol at S1 [ln ug/l]"
              ),
              missing = "no",
              missing_text = "Missing"
  )
#add p-value 
CAR.Tbl1 <- add_p(CAR.Tbl1,
                  test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "chisq.test"),
                  pvalue_fun = purrr::partial(style_pvalue, digits = 3)) %>% bold_p() %>% add_overall()


### Mid Preg --------------------------------------------------------
demo <- as.data.frame(cort_phenos[cort_phenos$pregstage == "II" & cort_phenos$CAR_trend %in% c("positive", 
                                                                                              "negative"),]) %>% 
  select(c("CAR_trend", 
           "CAR_iAUC", 
           "gestage_weeks",
           "waking_time", 
           "lg_waking_cort_concentration_ug_l",
           "peak_concentration",
           "peak_timeslot",
           "season",
           "caseVScontrol",
           "Maternal_Age_Years",
           "Parity",
           "Maternal_Education",
           "Maternal_Smoking_During_Pregnancy",
           "Maternal_Mental_Disorders",
           "Maternal_Diabetes_Disorders_anyVSnone",
           "Maternal_Hypertensive_Disorders_anyVSnone",
           "Maternal_Body_Mass_Index_in_Early_Pregnancy",
           "weight_gain_during_preg",
           "relative_to_recomd_weight_gain"))

CAR.Tbl2 <- 
  tbl_summary(data = demo, 
              by = CAR_trend,
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"), 
              digits = list(all_continuous() ~ c(2, 2, 2)),
              missing = "no",
              missing_text = "Missing"
  )
#add p-value 
CAR.Tbl2 <- add_p(CAR.Tbl2,
                  test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "chisq.test"),
                  pvalue_fun = purrr::partial(style_pvalue, digits = 3)) %>% bold_p() %>% add_overall()
#CAR.Tbl2 <- as_flex_table(CAR.Tbl2)
#save(CAR.Tbl2, file = "/Users/alicia_schowe/Desktop/ITU_final/Cortisol_Phenotypes/results/CAR.Tbl2.Rdata")

# gtsave(CAR.Tbl2, file = "Cortisol_Phenotypes/results/CAR_pregII.html")
# webshot(
#   url  = "Cortisol_Phenotypes/results/CAR_pregII.html", 
#   file = "Cortisol_Phenotypes/results/CAR_pregII.png",
#   zoom = 2   # doubles the resolution
#)
### Late Preg --------------------------------------------------------
demo <- as.data.frame(cort_phenos[cort_phenos$pregstage == "III" & cort_phenos$CAR_trend %in% c("positive", 
                                                                                               "negative"),]) %>% 
  select(c("CAR_trend", 
           "CAR_iAUC", 
           "gestage_weeks",
           "waking_time", 
           "lg_waking_cort_concentration_ug_l",
           "peak_concentration",
           "peak_timeslot",
           "season",
           "caseVScontrol",
           "Maternal_Age_Years",
           "Parity",
           "Maternal_Education",
           "Maternal_Smoking_During_Pregnancy",
           "Maternal_Mental_Disorders",
           "Maternal_Diabetes_Disorders_anyVSnone",
           "Maternal_Hypertensive_Disorders_anyVSnone",
           "Maternal_Body_Mass_Index_in_Early_Pregnancy",
           "weight_gain_during_preg",
           "relative_to_recomd_weight_gain"))

CAR.Tbl3 <- 
  tbl_summary(data = demo, 
              by = CAR_trend,
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"), 
              digits = list(all_continuous() ~ c(2, 2, 2)), 
              missing = "no",
              missing_text = "Missing"
  )
#add p-value 
CAR.Tbl3 <- add_p(CAR.Tbl3,
                  test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "chisq.test"),
                  pvalue_fun = purrr::partial(style_pvalue, digits = 3)) %>% bold_p() %>% add_overall()

remove(list=ls())

