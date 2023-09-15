## ALYS cohort AUCi 
# 24.08.2022

dat <- openxlsx::read.xlsx("Cortisol_Phenotypes/results/AYLS_AUCi.xlsx")
head(dat) #estimates provided by Polina Girchenko

# mean
mean(dat$AUCi, na.rm=T) # 1.179
sd(dat$AUCi, na.rm=T)

length(dat$AUCi[dat$AUCi<0 | dat$AUCi == 0])/length(dat$AUCi) #0.18%
1-0.179
length(dat$AUCi[dat$AUCi>0])/length(dat$AUCi)
                
       