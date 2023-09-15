#  Barcelona cohort CAR follow up 


# libraries ------------------------------
library(foreign)

# Data -----------------------------------
dat <- read.spss(file = "ITU_Cortisol_Publication_FEB2023/Data/22-11-22-Maternal-Cortisol-Database-Barcelona.sav", to.data.frame = TRUE)
head(dat)
summary(dat)

CAR_dat <- dat[,c("sample_id", 
                  "cortisol_1b_1t",
                  "time_cortisol_1b_1t",
                  "cortisol_2b_1t",
                  "time_cortisol_2b_1t",
                  "cortisol_1b_2t",
                  "time_cortisol_1b_2t",
                  "cortisol_2b_2t",
                  "time_cortisol_2b_2t",
                  "cortisol_1b_3t",
                  "time_cortisol_1b_3t",
                  "cortisol_2b_3t",
                  "time_cortisol_2b_3t",
                  "weeks_1t",
                  "weeks_2t",
                  "weeks_cortisol_3t")]


## early pregnancy -------
CAR_datI <- subset(CAR_dat, !is.na(cortisol_2b_1t) & !is.na(cortisol_1b_1t))
CAR_datI$diff_t1 <- CAR_datI$cortisol_2b_1t - CAR_datI$cortisol_1b_1t
summary(CAR_datI$diff_t1) #113
sd(CAR_datI$diff_t1)
length(CAR_datI$diff_t1[CAR_datI$diff_t1 <= 0 & !is.na(CAR_datI$diff_t1)])/length(CAR_datI$diff_t1[!is.na(CAR_datI$diff_t1)]) 
mean(CAR_datI$weeks_1t, na.rm=T)

## mid pregnancy -------
CAR_datII <- subset(CAR_dat, !is.na(cortisol_2b_2t) & !is.na(cortisol_1b_2t))
CAR_datII$diff_t2 <- CAR_datII$cortisol_2b_2t - CAR_datII$cortisol_1b_2t
nrow(CAR_datII)
summary(CAR_datII$diff_t2) #117
sd(CAR_datII$diff_t2)
length(CAR_datII$diff_t2[CAR_datII$diff_t2 <= 0 & !is.na(CAR_datII$diff_t2)])/length(CAR_datII$diff_t2[!is.na(CAR_datII$diff_t2)]) 

mean(CAR_datII$weeks_2t, na.rm=T)

## late pregnancy -------
CAR_datIII <- subset(CAR_dat, !is.na(cortisol_2b_3t) & !is.na(cortisol_1b_3t))
CAR_datIII$diff_t3 <- as.numeric(as.character(CAR_datIII$cortisol_2b_3t)) - as.numeric(as.character(CAR_datIII$cortisol_1b_3t))
summary(CAR_datIII$diff_t3) 
sd(CAR_datIII$diff_t3)
length(CAR_datIII$diff_t3[CAR_datIII$diff_t3 <= 0 & !is.na(CAR_datIII$diff_t3)])/length(CAR_datIII$diff_t3[!is.na(CAR_datIII$diff_t3)])

CAR_datIII$weeks_cortisol_3t <- as.numeric(as.character(CAR_datIII$weeks_cortisol_3t))
mean(as.numeric(CAR_datIII$weeks_cortisol_3t), na.rm=T)
