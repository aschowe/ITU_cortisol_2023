# ITU sample descriptives  ###############################################
# Purpose: Show sample characteristics
# 18.11.22 
# alicia_schowe@psych.mpg.de


# libraries -------------------------------------------------------------------
library(gtsummary)
library(dplyr)
library(gt)
library(naniar)
library(nlme)

# Data -------------------------------------------------------------------------
load("Data/ITU_combined_cortisol_dates_times_wide_format.Rdata")
load("Data/register_dat.final.Rdata")
characteristics <- subset(register_dat.final, participantID %in% wide_cort$participantID)

characteristics$Nulliparous <- characteristics$Parity
levels(characteristics$Nulliparous) <- c("yes","no")
characteristics$Chromosomal_Testing <- characteristics$caseVScontrol
levels(characteristics$Chromosomal_Testing) <- c("yes","no")
levels(characteristics$Maternal_Hypertensive_Disorders_anyVSnone) <- c("no","yes")
levels(characteristics$Maternal_Diabetes_Disorders_anyVSnone) <- c("no","yes")

# Analyses #############################################################

## Table 1: Sample Characteristics -------------------------------------- 
demo <- as.data.frame(characteristics) %>% select(
  colnames(characteristics[,c("Chromosomal_Testing",
                              "Maternal_Age_Years",
                              "Maternal_Education",
                              "Nulliparous",
                              "Maternal_Smoking_During_Pregnancy",
                              "Maternal_Mental_Disorders",
                              "Maternal_Body_Mass_Index_in_Early_Pregnancy",
                              "Maternal_Metabolic_Conditions",
                              "Maternal_Gestational_Diabetes",
                              "Maternal_Hypertensive_Disorders_anyVSnone" #,
                              )]) # n = 667 
)


Chromosomal_Testing.Tbl1 <- 
  tbl_summary(data = demo, 
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"), 
              digits = list(all_continuous() ~ c(2, 2, 2)),
              missing = "ifany",
              missing_text = "Missing"
              )




## Table 2: Cortisol Characteristics --------------------------------------------------
cort <- left_join(wide_cort, characteristics)
cort$CAR_absent <- ifelse(cort$CAR_iAUC <= 0, "yes", "no")
cort <- as.data.frame(cort)[,c("pregstage",
                       "gestage_weeks",
                       "waking_time",
                       "cort_concentration_ug_l_S1",
                       "cort_concentration_ug_l_S2",
                       "cort_concentration_ug_l_S3",
                       "cort_concentration_ug_l_S4",
                       "cort_concentration_ug_l_S5",
                       "cort_concentration_ug_l_S6",
                       "cort_concentration_ug_l_S7",
                       "CAR_iAUC",
                       "CAR_absent",
                       "cort_d_AUC")]

demo <- cort %>% select(
  colnames(cort)
)

Paper.Table.2 <- 
  tbl_summary(data = demo, 
              by = pregstage,
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"), 
              digits = list(all_continuous() ~ c(2, 2, 2)),           
              label = list(gestage_weeks ~ "Gestational Week",
                           CAR_iAUC ~ "CAR [AUCi]",
                           cort_d_AUC ~ "Diurnal Cortisol [AUCg]",
                           waking_time ~ "Time at Awakening",
                           CAR_absent ~ "Women with absent CAR",
                           cort_concentration_ug_l_S1 ~ "Cortisol Concentration S1 [ug/l]",
                           cort_concentration_ug_l_S2 ~ "Cortisol Concentration S2 [ug/l]",
                           cort_concentration_ug_l_S3 ~ "Cortisol Concentration S3 [ug/l]",
                           cort_concentration_ug_l_S4 ~ "Cortisol Concentration S4 [ug/l]",
                           cort_concentration_ug_l_S5 ~ "Cortisol Concentration S5 [ug/l]",
                           cort_concentration_ug_l_S6 ~ "Cortisol Concentration S6 [ug/l]",
                           cort_concentration_ug_l_S7 ~ "Cortisol Concentration S7 [ug/l]"),
              missing = "no",
              missing_text = "Missing"
              )



## adding ICC descriptives to Table 2 -----------------------------------------------------
cortisol_samples <- c( "waking_time",
                       "cort_concentration_ug_l_S1",
                       "cort_concentration_ug_l_S2",
                       "cort_concentration_ug_l_S3",
                       "cort_concentration_ug_l_S4",
                       "cort_concentration_ug_l_S5",
                       "cort_concentration_ug_l_S6",
                       "cort_concentration_ug_l_S7",
                       ) 

for(i in 1:length(cortisol_samples)){
  cort_index <- cortisol_samples[[i]]
S1 <- wide_cort[,c("participantID", "pregstage", cort_index)] 
S1 <- tidyr::spread(S1, pregstage, cort_index)
S1 <- S1[complete.cases(S1),]
#names(S1) <- paste0(names(S1), cort_index)
print(cor(S1[,c("I", "II", "III")]))
#print(cor.test(S1$II, S1$III))
print(irr::icc(S1[,c(2:4)], model = "twoway",
          type = "consistency", unit = "single"))
}
wide_cort$peak <- NA
for(i in 1:nrow(wide_cort)){
  wide_cort$peak[[i]] <- max(wide_cort$cort_concentration_ug_l_S1[[i]],
                             wide_cort$cort_concentration_ug_l_S2[[i]],
                             wide_cort$cort_concentration_ug_l_S3[[i]])
}
S1 <- wide_cort[,c("participantID", "pregstage", "peak")] 
S1 <- tidyr::spread(S1, pregstage, peak)
S1 <- S1[complete.cases(S1),]
print(irr::icc(S1[,c(2:4)], model = "twoway",
               type = "consistency", unit = "single"))

wide_cort %>% 
  group_by(pregstage) %>% 
  summarise(p=mean(peak,na.rm=T))

## ICC
AUCg <- wide_cort[,c("participantID", "pregstage", "cort_d_AUC")]
AUCg <- tidyr::spread(AUCg, pregstage, cort_d_AUC)
AUCg <- AUCg[complete.cases(AUCg),]
irr::icc(AUCg[,c(2:4)], model = "twoway",
    type = "consistency", unit = "single")


source("/Users/alicia_schowe/Desktop/ITU_final/Cortisol_Phenotypes/ITU_lmm_helper_functions.R")
ICC <- function(lme_model){
  icc <- as.numeric(VarCorr(lme_model)[4,1])/sum(as.numeric(VarCorr(lme_model)[c(2,4,5),1]))
  return(round(icc, digits = 2))
}
M.AUCg <- lme(cort_d_AUC ~
                waking_time +
                pregstage +
                season,
              random = list(~1|ELISA_analysis_plate,
                            ~ 1|participantID),
              method = "ML",
              data = wide_cort[wide_cort$participantID %in% AUCg$participantID,],
              na.action = "na.omit",
              control= lmeControl(msMaxIter = 200))
ICC(M.AUCg)

## ICC
AUCi <- wide_cort[,c("participantID", "pregstage", "CAR_iAUC")]
AUCi <- tidyr::spread(AUCi, pregstage, CAR_iAUC)
AUCi <- AUCi[complete.cases(AUCi),]
irr::icc(AUCi[,c(2:4)], model = "twoway",
         type = "consistency", unit = "single")

M.AUCi <- lme(CAR_iAUC ~
                waking_time +
                pregstage +
                season,
              random = list(~1|ELISA_analysis_plate,
                            ~ 1|participantID),
              method = "ML",
              data = wide_cort[wide_cort$participantID %in% AUCi$participantID,],
              na.action = "na.omit",
              control= lmeControl(msMaxIter = 200))
ICC(M.AUCi)

### assessment days 
wide_cort <- wide_cort %>% 
  group_by(participantID) %>%
  mutate(Nr_Stages = n())

Stages <- subset(wide_cort[,c("participantID", 
                     "Nr_Stages")], !duplicated(participantID))
summary(factor(Stages$Nr_Stages))

## table notes: absent CAR 

absent <- subset(wide_cort, CAR_trend == "negative" & Nr_Stages >= 2) 
length(unique(absent$participantID))

always_absent <- subset(wide_cort, CAR_trend_across_preg == "negative" & Nr_Stages >= 2)
length(unique(always_absent$participantID))

### document available number of  samples per women/pregnancy stage
#### Across gestation 
load("ITU_Cortisol_Publication_FEB2023/Data/ITU_combined_cortisol_dates_times_long_format.Rdata")
Samples <- long_cort %>% 
              group_by(participantID) %>% 
              mutate(nr_samples_overall = length(lg_cort_concentration_ug_l[!is.na(lg_cort_concentration_ug_l)]))

Samples <- subset(Samples[,c("participantID", 
                              "nr_samples_overall")], !duplicated(participantID))

mean_nr_singles_overall <- round(mean(Samples$nr_samples_overall), digits = 2)
range_nr_singles_overall <- range(Samples$nr_samples_overall)

#### Per Pregnancy Stage
PregstageSamples <- long_cort %>% 
  group_by(participantID, pregstage) %>% 
  mutate(nr_samples = length(lg_cort_concentration_ug_l[!is.na(lg_cort_concentration_ug_l)]))

PregstageSamples <- subset(PregstageSamples[,c("participantID", 
                             "nr_samples", "pregstage")], !duplicated(PregstageSamples[,c("participantID", 
                                                                                          "nr_samples", "pregstage")]))

mean_nr_singlesI <- round(mean(PregstageSamples$nr_samples[PregstageSamples$pregstage == "I"]), digits = 2)
range_nr_singlesI <- range(PregstageSamples$nr_samples[PregstageSamples$pregstage == "I"])
  
mean_nr_singlesII <- round(mean(PregstageSamples$nr_samples[PregstageSamples$pregstage == "II"]), digits = 2)
range_nr_singlesII <- range(PregstageSamples$nr_samples[PregstageSamples$pregstage == "II"])

mean_nr_singlesIII <- round(mean(PregstageSamples$nr_samples[PregstageSamples$pregstage == "III"]), digits = 2)
range_nr_singlesIII <- range(PregstageSamples$nr_samples[PregstageSamples$pregstage == "III"])

#### Complete Profiles 
percent_all7 <- round((length(unique(PregstageSamples$participantID[PregstageSamples$nr_samples == 7]))/N) *100, digits = 1)

# Supplemental Table 2: self-reported sampling times ========================
cort <- left_join(wide_cort, characteristics)
cort$CAR_absent <- ifelse(cort$CAR_iAUC <= 0, "yes", "no")
cort <- as.data.frame(cort)[,c("pregstage",
                               "time_of_day_in_hours_since_midnight_S1",
                               "time_of_day_in_hours_since_midnight_S2",
                               "time_of_day_in_hours_since_midnight_S3",
                               "time_of_day_in_hours_since_midnight_S4",
                               "time_of_day_in_hours_since_midnight_S5",
                               "time_of_day_in_hours_since_midnight_S6",
                               "time_of_day_in_hours_since_midnight_S7")]

demo <- cort %>% select(
  colnames(cort)
)

Paper.Table.2b <- 
  tbl_summary(data = demo, 
              by = pregstage,
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"), 
              digits = list(all_continuous() ~ c(2, 2, 2)),           
              label = list(
              )
  )



## Supplementary Sample Comparisons:  All women vs women without cortisol data -------------------
library(naniar)
#### Load data 
load("/Users/alicia_schowe/Desktop/ITU_final/ITU_Cortisol_Publication_FEB2023/Data/register_dat.final.Rdata")

maternal_edu <- read.delim("/Users/alicia_schowe/Desktop/ITU_final/ITU_Cortisol_Publication_FEB2023/Data/ITU maternal education.dat")
names(maternal_edu)[1] <- "participantID"
names(maternal_edu)[3] <- "Maternal_Education"
maternal_edu <- maternal_edu[,c("participantID",
                                "Maternal_Education")]

maternal_edu <- maternal_edu %>% replace_with_na(replace = list(Maternal_Education = -9))

maternal_edu$Maternal_Education <- factor(maternal_edu$Maternal_Education,
                                          levels = c(1,2,3),
                                          labels = c("primary", "applied university", "university"))

medication <- read.delim("/Users/alicia_schowe/Desktop/ITU_final/Cortisol_Processing/Data/ITU_psychotrophicmedication05July21_Maternal_CurrentPregnancy.dat")
names(medication)[1] <- "participantID"


#### Merge all 
Complete_ITU <- subset(register_dat.final, participantID %in% medication$participantID)
Complete_ITU <- left_join(Complete_ITU,
                          maternal_edu)

Complete_ITU$cort_avail <- NA
Complete_ITU$cort_avail[Complete_ITU$participantID %in% wide_cort$participantID] <- "yes"
Complete_ITU$cort_avail[!Complete_ITU$participantID %in% wide_cort$participantID] <- "no"

#### Compare all by cortisol availability  

demo <- as.data.frame(Complete_ITU) %>% select(
  colnames(Complete_ITU[,-c(1)]) # n = 955
)

demo$cort_avail <- factor(demo$cort_avail)

Supp_Table.cortAvail <- 
  tbl_summary(data = demo, 
              by = cort_avail,
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"), 
              digits = list(all_continuous() ~ c(2, 2, 2)),           
              label = list(
  ))
#add p-value 
Supp_Table.cortAvail <- add_p(Supp_Table.cortAvail,
                              test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "chisq.test"),
                              pvalue_fun = purrr::partial(style_pvalue, digits = 3)) %>% bold_p() %>% add_overall()
