# ITU: Plots 
# 18.11.22
# alicia_schowe@psych.mpg.de

library(ggplot2)
library(dplyr)
library(patchwork)


# Data -----------------------------------------------------------
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
cort_phenos_wide <- left_join(wide_cort, register_dat.final)
levels(cort_phenos_wide$caseVScontrol) <- c("Yes", "No")

load("Data/CAR.df.Rdata")

load("Data/diurnal.df.Rdata")
levels(diurnal.df$pregstage) <- c("Early", "Mid", "Late")

for(i in 1:nrow(diurnal.df)){
  S <- diurnal.df$timeslot[[i]]
  if(S %in% c("S1", "S2", "S3")){
    diurnal.df$timeslot[[i]] <- "peak"
  }
  else{diurnal.df$timeslot[[i]] <- S}
}

# Figure 1 -------------------------------------------------------

DIUR <- diurnal.df %>% 
  group_by(pregstage, timeslot) %>% 
  summarise(
    n = n(),
    c = mean(lg_cort_concentration_ug_l, na.rm = T),
    SE = sd(lg_cort_concentration_ug_l, na.rm = T)/sqrt(n),
    t = mean(time_of_day_in_hours_since_midnight, na.rm = T))

theme_set(theme_classic())
DIUR <- ggplot(na.omit(DIUR), 
       aes(x = t,
           y = c,
           colour = pregstage,
           shape = pregstage,
           linetype = pregstage)) +
  geom_line(aes(group = pregstage)) + 
  geom_line(data = na.omit(DIUR)) + 
  geom_point(aes(x = t, shape = pregstage), size = 2) +
  geom_errorbar(aes(ymin=c-SE, ymax=c+SE), width=.1) +
  ylab("Cortisol Concentration [Ln(ug/l+1)]") +
  scale_x_continuous(name="", 
                     breaks = c(7.5, 12, 17, 22.5),
                     labels = c("7:30", "12:00", "17:00", "22:30")) + 
  theme(text = element_text(size=12)) + 
  theme(legend.position="bottom",) +
  labs(color='Pregnancy Stage', linetype = "Pregnancy Stage", shape = "Pregnancy Stage") 
  

cort_phenos_wide$pregstage <- factor(cort_phenos_wide$pregstage)
levels(cort_phenos_wide$pregstage) <- c("Early", "Mid", "Late")

bxp_AUCg <- ggpubr::ggboxplot(
  cort_phenos_wide, x = "pregstage", y = "cort_d_AUC",
  color = "pregstage",
  xlab = "",  
  #title = "Cortisol Concentration by Pregnancy Stage,
  legend = "none",
  legend.title = NULL,
  ylab = "Cortsiol Concentration [AUCg]"
)

setEPS()
postscript("Schowe_2023_fig1.eps")
DIUR + inset_element(bxp_AUCg, left = 0.5, bottom = 0.5, right = 1, top = 1)
dev.off()

# ============================================================================
# Figure 2: CAR and BMI -------------------------------------------------------------
CAR.df$BMI_highvslow <- ifelse(CAR.df$Maternal_Body_Mass_Index_in_Early_Pregnancy >= 25, 
                               "BMI >= 25", 
                               "BMI < 25")
CAR.df$BMI_highvslow <- factor(CAR.df$BMI_highvslow,
                               ordered = FALSE)

CAR.df$BMI_highvslow <- relevel(CAR.df$BMI_highvslow,
                                ref = "BMI >= 25")

levels(CAR.df$pregstage) <- c("Early", "Mid", "Late")

CAR <- CAR.df %>% 
  group_by(timeslot, BMI_highvslow) %>% 
  summarise(
    n = n(),
    c = mean(lg_cort_concentration_ug_l, na.rm = T),
    SE = sd(lg_cort_concentration_ug_l, na.rm = T)/sqrt(n),
    t = mean(hours_since_waking, na.rm = T))

theme_set(theme_classic())
P <- ggplot(na.omit(CAR), 
            aes(x = t,
                y = c,
                colour = BMI_highvslow,
                linetype = BMI_highvslow,
                shape = BMI_highvslow)) +
  geom_line(aes(group = BMI_highvslow)) + 
  #facet_grid(. ~ pregstage) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5)) + 
  geom_line(data = na.omit(CAR)) + 
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=c-SE, ymax=c+SE), width=.1) +
  ylab("Cortisol Concentration [Ln(ug/l+1)]") +
  scale_x_continuous(name="", 
                     breaks = c(0,0.25,0.5),
                     labels = c("Awakening", "+15", "+30")) + 
  theme(text = element_text(size=12),
        axis.title.x=element_blank()) + 
  labs(color='Body Mass Index', linetype = "Body Mass Index", shape = "Body Mass Index") + 
  theme(legend.position="bottom") 

P


setEPS()
postscript("Schowe_2023_fig2.eps")
P 
dev.off()

# =============================================================================
# Figure 3: DCS and BMI 
diurnal.df$BMI_highvslow <- ifelse(diurnal.df$Maternal_Body_Mass_Index_in_Early_Pregnancy >= 25, 
                               "BMI >= 25", 
                               "BMI < 25")
diurnal.df$BMI_highvslow <- factor(diurnal.df$BMI_highvslow,
                               ordered = FALSE)

diurnal.df$BMI_highvslow <- relevel(diurnal.df$BMI_highvslow,
                                ref = "BMI < 25")

DIUR <- diurnal.df %>% 
  group_by(pregstage, timeslot, BMI_highvslow) %>% 
  summarise(
    n = n(),
    c = mean(lg_cort_concentration_ug_l, na.rm = T),
    SE = sd(lg_cort_concentration_ug_l, na.rm = T)/sqrt(n),
    t = mean(time_of_day_in_hours_since_midnight, na.rm = T))

theme_set(theme_classic())
DIUR <- ggplot(na.omit(DIUR), 
               aes(x = t,
                   y = c,
                   colour = pregstage,
                   linetype = pregstage, 
                   shape = pregstage)) +
  geom_line(aes(group = pregstage)) + 
  facet_grid(. ~ BMI_highvslow) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5)) + 
  geom_line(data = na.omit(DIUR)) + 
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=c-SE, ymax=c+SE), width=.1) +
  ylab("Cortisol Concentration [Ln(ug/l+1)]") +
  scale_x_continuous(name="", 
                     breaks = c(7.5, 12, 17, 22.5),
                     labels = c("7:30", "12:00", "17:00", "22:30")) + 
  theme(text = element_text(size=12),
        axis.title.x=element_blank()) + 
  labs(color='Pregnancy Stage', linetype = "Pregnancy Stage", shape = "Pregnancy Stage") + 
  theme(legend.position="bottom") 

setEPS()
postscript("Schowe_2023_fig3.eps")
DIUR
dev.off()

