# ITU: Exploratory Analyses between maternal characteristics and CAR  ############
# 18.11.22
# alicia_schowe@pysch.mpg.de 

#setwd("ITU_Cortisol_Publication_FEB2023/")
# Data -----------------------------------------------------------------------
load("Data/CAR.df.Rdata")

# libraries -------------------------------------------------------------------
library(nlme)
library(data.table)
library(APAstyler) #remotes::install_github("darrellpenta/APAstyler")
library(naniar)
library(dplyr)

# Functions -------------------------------------------------------------------
source("ITU_lmm_helper_functions.R")

# Cardiometabolic variables of interest ---------------------------------------
phenos <- c(
  "Maternal_Hypertensive_Disorders_anyVSnone",
  "Maternal_Gestational_Diabetes",
  "Maternal_Body_Mass_Index_in_Early_Pregnancy_cent",
  "Maternal_Metabolic_Conditions"
  )

# CAR Basic Model --------------------------------------------------------------
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


# Separate effects -----------------------------------------------------------
CAR.results <- data.frame(matrix(ncol = 9, nrow = 0))
names(CAR.results) <- c("M1_b", "M1_SE", "M1_t", "M1_p", 
                        "M2_b", "M2_SE", "M2_t", "M2_p",
                        "M2_Interpretation")

covariates <- c("caseVScontrol",
                "Parity",
                "Maternal_Age_Years",
                "Maternal_Education",
                "Maternal_Smoking_During_Pregnancy",
                "Maternal_Mental_Disorders")

## Main Effect
for(i in 1:length(phenos)){
  pheno <- phenos[[i]]
  if(!pheno ==  "Maternal_Metabolic_Conditions"){
    M1_form <- paste0(".~.", "+", pheno) 
    M1 <- update(M.CAR.full, M1_form)
    RES.CAR_final_M1 <- model_res(summary(M1))
    RES.CAR_final_M1 <- mutate_all(RES.CAR_final_M1,funs(replace(., is.na(.), "")))
    CAR.results[pheno,c(1:4)] <- RES.CAR_final_M1[grep(pheno, rownames(RES.CAR_final_M1)), c(1:4)]
  
    #add covariats
    M2_form <- paste0(".~.", "+", pheno, "+", paste(covariates, collapse = "+")) 
    M2 <- update(M.CAR.full, M2_form)
    RES.CAR_final_M2 <- model_res(summary(M2))
    RES.CAR_final_M2 <- mutate_all(RES.CAR_final_M2,funs(replace(., is.na(.), "")))
    CAR.results[pheno,c(5:9)] <- RES.CAR_final_M2[grep(pheno, rownames(RES.CAR_final_M2)), c(1:5)]
  }
  if(pheno == "Maternal_Metabolic_Conditions"){
    M1_form <- paste0(".~.", "+", pheno) 
    M1 <- update(M.CAR.full, M1_form)
    RES.CAR_final_M1 <- model_res(summary(M1))
    RES.CAR_final_M1 <- mutate_all(RES.CAR_final_M1,funs(replace(., is.na(.), "")))
    #print(RES.CAR_final_M1)
    CAR.results["Maternal_Metabolic_Conditions1",c(1:4)] <- 
      RES.CAR_final_M1[grep("Maternal_Metabolic_Conditions1", rownames(RES.CAR_final_M1)), c(1:4)]
    CAR.results["Maternal_Metabolic_Conditions>1",c(1:4)] <- 
      RES.CAR_final_M1[grep("Maternal_Metabolic_Conditions>1", rownames(RES.CAR_final_M1)), c(1:4)]
    
    #add covariats
    M2_form <- paste0(".~.", "+", pheno, "+", paste(covariates, collapse = "+"))
    M2 <- update(M.CAR.full, M2_form)
    RES.CAR_final_M2 <- model_res(summary(M2))
    RES.CAR_final_M2 <- mutate_all(RES.CAR_final_M2,funs(replace(., is.na(.), "")))
    CAR.results["Maternal_Metabolic_Conditions1",c(5:9)] <- 
      RES.CAR_final_M2[grep("Maternal_Metabolic_Conditions1", rownames(RES.CAR_final_M2)), c(1:5)]
    CAR.results["Maternal_Metabolic_Conditions>1",c(5:9)] <- 
      RES.CAR_final_M2[grep("Maternal_Metabolic_Conditions>1", rownames(RES.CAR_final_M2)), c(1:5)]
  }
}

## by Slope
for(i in 1:length(phenos)){
  pheno <- phenos[[i]]
    M1_form <- paste0(".~.", "+", pheno, "+", pheno, ": hours_since_waking") 
    M1 <- update(M.CAR.full, M1_form)
    RES.CAR_final_M1 <- model_res(summary(M1))
    RES.CAR_final_M1 <- mutate_all(RES.CAR_final_M1,funs(replace(., is.na(.), "")))
    if(!pheno ==  "Maternal_Metabolic_Conditions"){
      label <- paste0(pheno, "x Time")
      CAR.results[label,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":", pheno), rownames(RES.CAR_final_M1)), c(1:4)]
    } else{
        label1 <- paste0("Maternal_Metabolic_Conditions1", "x Time")
        label2 <- paste0("Maternal_Metabolic_Conditions>1", "x Time")
        CAR.results[label1,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M1)), c(1:4)]
        CAR.results[label2,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M1)), c(1:4)]
    }
    #add covariats
    M2_form <- paste0(".~.", "+", pheno, "+", 
                     "+", pheno, ": hours_since_waking", "+",
                     paste(covariates, collapse = "+")) 
    M2 <- update(M.CAR.full, M2_form)
    RES.CAR_final_M2 <- model_res(summary(M2))
    RES.CAR_final_M2 <- mutate_all(RES.CAR_final_M2,funs(replace(., is.na(.), "")))
    if(!pheno ==  "Maternal_Metabolic_Conditions"){
    CAR.results[label,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":", pheno), rownames(RES.CAR_final_M2)), c(1:5)]
    } else{ 
      CAR.results[label1,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M2)), c(1:5)]
      CAR.results[label2,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M2)), c(1:5)]
    }
}

## effect dependent on pregnancy stage 

### pregnancy stage * Concentration 
for(i in 1:length(phenos)){
  pheno <- phenos[[i]]
    M1_form <- paste0(".~.", "+", pheno, "+", 
                      pheno, ": pregstage") 
    M1 <- update(M.CAR.full, M1_form)
   RES.CAR_final_M1 <- model_res(summary(M1))
    RES.CAR_final_M1 <- mutate_all(RES.CAR_final_M1,funs(replace(., is.na(.), "")))
    if(!pheno ==  "Maternal_Metabolic_Conditions"){
      labelT2 <- paste(pheno, "x T2 [ref = T1]")
      labelT3 <- paste(pheno, "x T3 [ref = T1]")
  
      CAR.results[labelT2,c(1:4)] <- RES.CAR_final_M1[grep(paste0("T2:", pheno), rownames(RES.CAR_final_M1)), c(1:4)]
      CAR.results[labelT3,c(1:4)] <- RES.CAR_final_M1[grep(paste0("T3:", pheno), rownames(RES.CAR_final_M1)), c(1:4)]
    }else{
      labelT2a <- paste("Maternal_Metabolic_Conditions1 x T2 [ref = T1]")
      labelT3a <- paste("Maternal_Metabolic_Conditions1 x T3 [ref = T1]")
      CAR.results[labelT2a,c(1:4)] <- RES.CAR_final_M1[grep(paste0("T2:", "Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M1)), c(1:4)]
      CAR.results[labelT3a,c(1:4)] <- RES.CAR_final_M1[grep(paste0("T3:", "Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M1)), c(1:4)]
      
      labelT2b <- paste("Maternal_Metabolic_Conditions>1 x T2 [ref = T1]")
      labelT3b <- paste("Maternal_Metabolic_Conditions>1 x T3 [ref = T1]")
      CAR.results[labelT2b,c(1:4)] <- RES.CAR_final_M1[grep(paste0("T2:", "Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M1)), c(1:4)]
      CAR.results[labelT3b,c(1:4)] <- RES.CAR_final_M1[grep(paste0("T3:", "Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M1)), c(1:4)]
    }
    #add covariats
    M2_form <- paste0(".~.", "+", pheno, "+", 
                      "+", pheno, ": pregstage", "+",
                     paste(covariates, collapse = "+")) 
   M2 <- update(M.CAR.full, M2_form)
   RES.CAR_final_M2 <- model_res(summary(M2))
   RES.CAR_final_M2 <- mutate_all(RES.CAR_final_M2,funs(replace(., is.na(.), "")))
   if(!pheno ==  "Maternal_Metabolic_Conditions"){
    CAR.results[labelT2,c(5:9)] <- RES.CAR_final_M2[grep(paste0("T2:", pheno), rownames(RES.CAR_final_M2)), c(1:5)]
    CAR.results[labelT3,c(5:9)] <- RES.CAR_final_M2[grep(paste0("T3:", pheno), rownames(RES.CAR_final_M2)), c(1:5)]
   }else{
     CAR.results[labelT2a,c(5:9)] <- RES.CAR_final_M2[grep(paste0("T2:", "Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M2)), c(1:5)]
     CAR.results[labelT3a,c(5:9)] <- RES.CAR_final_M2[grep(paste0("T3:", "Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M2)), c(1:5)]
     CAR.results[labelT2b,c(5:9)] <- RES.CAR_final_M2[grep(paste0("T2:", "Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M2)), c(1:5)]
     CAR.results[labelT3b,c(5:9)] <- RES.CAR_final_M2[grep(paste0("T3:", "Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M2)), c(1:5)]
   }
}

### pregnancy stage * Time
for(i in 1:length(phenos)){
  pheno <- phenos[[i]]
    M1_form <- paste0(".~.", "+", pheno, "+", 
                      pheno, "*hours_since_waking*pregstage") 
    M1 <- update(M.CAR.full, M1_form)
    RES.CAR_final_M1 <- model_res(summary(M1))
    RES.CAR_final_M1 <- mutate_all(RES.CAR_final_M1,funs(replace(., is.na(.), "")))
    
    if(!pheno ==  "Maternal_Metabolic_Conditions"){
    labelT2 <- paste(pheno, "x Time x T2 [ref = T1]")
    labelT3 <- paste(pheno, "x Time x T3 [ref = T1]")
    CAR.results[labelT2,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":pregstageT2:", pheno), rownames(RES.CAR_final_M1)), c(1:4)]
    CAR.results[labelT3,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":pregstageT3:", pheno), rownames(RES.CAR_final_M1)), c(1:4)]  
    }else{
      labelT2a <- paste("Maternal_Metabolic_Conditions1 x Time x T2 [ref = T1]")
      labelT3a <- paste("Maternal_Metabolic_Conditions1", "x Time x T3 [ref = T1]")
      CAR.results[labelT2a,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":pregstageT2:Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M1)), c(1:4)]
      CAR.results[labelT3a,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":pregstageT3:Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M1)), c(1:4)]
      
      labelT2b <- paste("Maternal_Metabolic_Conditions>1 x Time x T2 [ref = T1]")
      labelT3b <- paste("Maternal_Metabolic_Conditions>1", "x Time x T3 [ref = T1]")
      CAR.results[labelT2b,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":pregstageT2:Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M1)), c(1:4)]
      CAR.results[labelT3b,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":pregstageT3:Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M1)), c(1:4)]
    }
    #add covariats
    M2_form <- paste0(".~.", "+", pheno, "+",
                    pheno, "*hours_since_waking*pregstage", "+",
                    paste(covariates, collapse = "+"))
    M2 <- update(M.CAR.full, M2_form)
    RES.CAR_final_M2 <- model_res(summary(M2))
    RES.CAR_final_M2 <- mutate_all(RES.CAR_final_M2,funs(replace(., is.na(.), "")))
    if(!pheno == "Maternal_Metabolic_Conditions"){
    CAR.results[labelT2,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":pregstageT2:", pheno), rownames(RES.CAR_final_M2)), c(1:5)]
    CAR.results[labelT3,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":pregstageT3:", pheno), rownames(RES.CAR_final_M2)), c(1:5)]
    }else{
      CAR.results[labelT2a,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":pregstageT2:Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M2)), c(1:5)]
      CAR.results[labelT3a,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":pregstageT3:Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M2)), c(1:5)]
      CAR.results[labelT2b,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":pregstageT2:Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M2)), c(1:5)]
      CAR.results[labelT3b,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":pregstageT3:Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M2)), c(1:5)]
  }
}



### save -----------------------------------------------------------------------
CAR.results <- tibble::rownames_to_column(CAR.results, "Predictor")

#save(CAR.results, file = "results/Tables/CAR.results_02052023.Rdata")
writexl::write_xlsx(CAR.results, "results/Tables/CAR.results_0252023.xlsx")

# relevel for supplementary tables --------------------------------------------
CAR.df$pregstage <- relevel(CAR.df$pregstage, ref = "T2")

CAR.resultsT2 <- data.frame(matrix(ncol = 9, nrow = 0))
names(CAR.resultsT2) <- c("M1_b", "M1_SE", "M1_t", "M1_p", 
                        "M2_b", "M2_SE", "M2_t", "M2_p",
                        "M2_Interpretation")

### pregnancy stage * Concentration 
for(i in 1:length(phenos)){
  pheno <- phenos[[i]]

    M1_form <- paste0(".~.", "+", pheno, "+", 
                      pheno, ": pregstage") 
    M1 <- update(M.CAR.full, M1_form)
    RES.CAR_final_M1 <- model_res(summary(M1))
    RES.CAR_final_M1 <- mutate_all(RES.CAR_final_M1,funs(replace(., is.na(.), "")))
    if(!pheno ==  "Maternal_Metabolic_Conditions"){
      labelT2 <- paste(pheno, "x T1 [ref = T2]")
      labelT3 <- paste(pheno, "x T3 [ref = T2]")
      
      CAR.resultsT2[labelT2,c(1:4)] <- RES.CAR_final_M1[grep(paste0("T1:", pheno), rownames(RES.CAR_final_M1)), c(1:4)]
      CAR.resultsT2[labelT3,c(1:4)] <- RES.CAR_final_M1[grep(paste0("T3:", pheno), rownames(RES.CAR_final_M1)), c(1:4)]
    }else{
      labelT2a <- paste("Maternal_Metabolic_Conditions1 x T1 [ref = T2]")
      labelT3a <- paste("Maternal_Metabolic_Conditions1 x T3 [ref = T2]")
      CAR.resultsT2[labelT2a,c(1:4)] <- RES.CAR_final_M1[grep(paste0("T1:", "Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M1)), c(1:4)]
      CAR.resultsT2[labelT3a,c(1:4)] <- RES.CAR_final_M1[grep(paste0("T3:", "Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M1)), c(1:4)]
      
      labelT2b <- paste("Maternal_Metabolic_Conditions>1 x T1 [ref = T2]")
      labelT3b <- paste("Maternal_Metabolic_Conditions>1 x T3 [ref = T2]")
      CAR.resultsT2[labelT2b,c(1:4)] <- RES.CAR_final_M1[grep(paste0("T1:", "Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M1)), c(1:4)]
      CAR.resultsT2[labelT3b,c(1:4)] <- RES.CAR_final_M1[grep(paste0("T3:", "Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M1)), c(1:4)]
    }
    #add covariates
    M2_form <- paste0(".~.", "+", pheno, "+", 
                      "+", pheno, ": pregstage", "+",
                      paste(covariates, collapse = "+")) 
    M2 <- update(M.CAR.full, M2_form)
    RES.CAR_final_M2 <- model_res(summary(M2))
    RES.CAR_final_M2 <- mutate_all(RES.CAR_final_M2,funs(replace(., is.na(.), "")))
    if(!pheno ==  "Maternal_Metabolic_Conditions"){
      CAR.resultsT2[labelT2,c(5:9)] <- RES.CAR_final_M2[grep(paste0("T1:", pheno), rownames(RES.CAR_final_M2)), c(1:5)]
      CAR.resultsT2[labelT3,c(5:9)] <- RES.CAR_final_M2[grep(paste0("T3:", pheno), rownames(RES.CAR_final_M2)), c(1:5)]
    }else{
      CAR.resultsT2[labelT2a,c(5:9)] <- RES.CAR_final_M2[grep(paste0("T1:", "Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M2)), c(1:5)]
      CAR.resultsT2[labelT3a,c(5:9)] <- RES.CAR_final_M2[grep(paste0("T3:", "Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M2)), c(1:5)]
      CAR.resultsT2[labelT2b,c(5:9)] <- RES.CAR_final_M2[grep(paste0("T1:", "Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M2)), c(1:5)]
      CAR.resultsT2[labelT3b,c(5:9)] <- RES.CAR_final_M2[grep(paste0("T3:", "Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M2)), c(1:5)]
    }
}

### pregnancy stage * Time
for(i in 1:length(phenos)){
  pheno <- phenos[[i]]
  M1_form <- paste0(".~.", "+", pheno, "+", 
                    pheno, "*hours_since_waking*pregstage") 
  M1 <- update(M.CAR.full, M1_form)
  RES.CAR_final_M1 <- model_res(summary(M1))
  RES.CAR_final_M1 <- mutate_all(RES.CAR_final_M1,funs(replace(., is.na(.), "")))
  
  if(!pheno ==  "Maternal_Metabolic_Conditions"){
    labelT2 <- paste(pheno, "x Time x T1 [ref = T2]")
    labelT3 <- paste(pheno, "x Time x T3 [ref = T2]")
    CAR.resultsT2[labelT2,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":pregstageT1:", pheno), rownames(RES.CAR_final_M1)), c(1:4)]
    CAR.resultsT2[labelT3,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":pregstageT3:", pheno), rownames(RES.CAR_final_M1)), c(1:4)]  
  }else{
    labelT2a <- paste("Maternal_Metabolic_Conditions1 x Time x T1 [ref = T2]")
    labelT3a <- paste("Maternal_Metabolic_Conditions1", "x Time x T3 [ref = T2]")
    CAR.resultsT2[labelT2a,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":pregstageT1:Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M1)), c(1:4)]
    CAR.resultsT2[labelT3a,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":pregstageT3:Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M1)), c(1:4)]
    
    labelT2b <- paste("Maternal_Metabolic_Conditions>1 x Time x T1 [ref = T2]")
    labelT3b <- paste("Maternal_Metabolic_Conditions>1", "x Time x T3 [ref = T2]")
    CAR.resultsT2[labelT2b,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":pregstageT1:Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M1)), c(1:4)]
    CAR.resultsT2[labelT3b,c(1:4)] <- RES.CAR_final_M1[grep(paste0(":pregstageT3:Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M1)), c(1:4)]
  }
  #add covariats
  M2_form <- paste0(".~.", "+", pheno, "+",
                    pheno, "*hours_since_waking*pregstage", "+",
                    paste(covariates, collapse = "+"))
  M2 <- update(M.CAR.full, M2_form)
  RES.CAR_final_M2 <- model_res(summary(M2))
  RES.CAR_final_M2 <- mutate_all(RES.CAR_final_M2,funs(replace(., is.na(.), "")))
  if(!pheno == "Maternal_Metabolic_Conditions"){
    CAR.resultsT2[labelT2,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":pregstageT1:", pheno), rownames(RES.CAR_final_M2)), c(1:5)]
    CAR.resultsT2[labelT3,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":pregstageT3:", pheno), rownames(RES.CAR_final_M2)), c(1:5)]
  }else{
    CAR.resultsT2[labelT2a,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":pregstageT1:Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M2)), c(1:5)]
    CAR.resultsT2[labelT3a,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":pregstageT3:Maternal_Metabolic_Conditions1"), rownames(RES.CAR_final_M2)), c(1:5)]
    CAR.resultsT2[labelT2b,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":pregstageT1:Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M2)), c(1:5)]
    CAR.resultsT2[labelT3b,c(5:9)] <- RES.CAR_final_M2[grep(paste0(":pregstageT3:Maternal_Metabolic_Conditions>1"), rownames(RES.CAR_final_M2)), c(1:5)]
  }
}
## save ---------------------------------------
CAR.resultsT2 <- tibble::rownames_to_column(CAR.resultsT2, "Predictor")

#save(CAR.resultsT2, file = "results/Tables/CAR.resultsT2_02052023.Rdata")
writexl::write_xlsx(CAR.resultsT2, "results/Tables/CAR.resultsT2_02052023.xlsx")
