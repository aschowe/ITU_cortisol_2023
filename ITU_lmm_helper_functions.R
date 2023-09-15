# Helper Functions to save lmm output 
# 01.03.22 
# alicia_schowe@gmx.de 

#libraries
library(APAstyler)

## Functions to compute ICC from 4-level mixed models 
### ICC
ICC_across <- function(lme_model){
  icc <- sum(as.numeric(VarCorr(lme_model)[c(4),1]))/sum(as.numeric(VarCorr(lme_model)[c(2,4,5,7,8),1]))
  return(round(icc, digits = 2))
}

ICC_within <- function(lme_model){
  icc <- sum(as.numeric(VarCorr(lme_model)[c(4,5,7),1]))/sum(as.numeric(VarCorr(lme_model)[c(2,4,5,7,8),1]))
  return(round(icc, digits = 2))
}

## Functions to save results ---------------------------------------------------

# round output of lmm tTable 
round_df <- function(x, cols1, digits1, cols2, digits2){
  x <- as.data.frame(x)
  # round all  variables
  # x: data frame 
  # cols: select which columns to round 
  # digits: number of digits to round corresponding to cols 
  for(i in 1:ncol(x[,c(cols1)])){
    colname <- names(x[,c(cols1)])[[i]]
    col <- x[,c(colname)]
    for(j in 1:length(col)){
      value <- as.numeric(col[[j]])
      if(!is.na(value)){
        # print(value)
        if(abs(value) > 0.005){
          x[,c(cols1)][[j,i]] <-  format(round(value, digits1))
        }
        else{x[,c(cols1)][[j,i]] <- format(value, digits=1, scientific = T)}
      }
    }
  }
  p_col <- x[,c(cols2)]
  for(k in 1:length(p_col)){
    value <- as.numeric(p_col[[k]])
    if(!is.na(value)){
      if(abs(value) >= 0.0001){
        x[k,c(cols2)] <-  format(round(value, digits2))
      }
      if(abs(value) >= 0.0001 & abs(value) < 0.0005){
        x[k,c(cols2)] <-  format(round(value, 4))
      }
      if(abs(value) < 0.0001){
        x[k,c(cols2)] <-  "<.0001"
      }
    }
  }
  return(x)
}

#calculate percent changes from beta coefficients 
res_Bchange <- function(res_df){ ##tTable of a summary model object as dataframe 
  res_df[["pchange"]] <- NA
  for(i in 1:nrow(res_df)){
    p <- as.numeric(res_df[["p-value"]][[i]])
    val <- res_df[["Value"]][[i]]
    if(p <= 0.05){
      res_df[["pchange"]][[i]] <- ((exp(val))-1)*100
    }
  }
  df <- res_df[,-c(3)]
  df[1,5] <- NA
  names(df) <- c("Coefficient", "SE", "t-statistic", "p", "Interpretation")
  df <- round_df(df, c(1:3,5), 3, 4,3)
  #to follow APA style, delete leading zero
  for(j in 1:nrow(df)){
    P <- df[[j, "p"]]
    if(P < 1 & P > 0){
      P <- snip(P, lead = 1)
      df[[j, "p"]] <- snip(P, lead = 1)
    }
  }
  return(df)
}

#extract the random part and compute R explained 
res_random <- function(model_summary){ #model summary object 
  random <- VarCorr(model_summary)[(c(2,4:5,7:8)), c(1,3)] #only works for 3-level model 
  random <- as.data.frame(random)
  names(random) <- c("Variance_Component", "Covariance")
  random$Variance_Component <- round(as.numeric(random$Variance_Component), digits = 3)
  random$Covariance <- as.numeric(random$Covariance)
  #random <- round(random, 1:2,2, 1:2,2)
  n <- nobs(model_summary)
  random <- rbind(random, c(n, NA))
  varests <- as.numeric(VarCorr(model_summary)[(c(2,4:5,7:8)), c(1)])
  total_R_explained <- format(round(sum(varests[1:4])/sum(varests)*100,1))
  random <- rbind(random, c(total_R_explained, NA))
  row.names(random) <- c("Intercept (ELISA_plate)",
                         "Intercept (person-level)", 
                         "Time", 
                         "Intercept (pregstage-level)", 
                         "Residual", 
                         "Number of Observations",
                         "Total R explained")
  return(random)
}

# extract fixed part and add percent explained 
model_res <- function(model_summary){
  tTable <- res_Bchange(as.data.frame(model_summary[["tTable"]]))
  return(tTable)
}

#extract fixed and random part together 
model_res.wRandom <- function(model_summary){
  tTable <- res_Bchange(as.data.frame(model_summary[["tTable"]]))
  random_part <- res_random(model_summary)
  random_part$empty1 <- NA
  random_part$empty2 <- NA
  random_part$empty3 <- NA
  tTable <- rbind(tTable,c("Variance Component", "Covariance", rep(NA,3)))
  names(random_part) <- names(tTable)
  Res <- rbind(tTable, random_part)
  Res <- setDT(Res, keep.rownames = "predictor")
  Res$predictor[rownames(Res) == Res$predictor] <- "Random Component"
  return(Res)
}

## formular for 3-level structure 
res_randomPeak <- function(model_summary){ #model summary object 
  random <- VarCorr(model_summary)[(c(2,4:5)), c(1:2)] #only works for 3-level model 
  random <- as.data.frame(random)
  names(random) <- c("Variance_Component", NA)
  random$Variance_Component <- round(as.numeric(random$Variance_Component), digits = 3)
  #random$Covariance <- as.numeric(random$Covariance)
  #random <- round(random, 1:2,2, 1:2,2)
  n <- nobs(model_summary)
  random <- rbind(random, c(n))
  varests <- as.numeric(VarCorr(model_summary)[(c(2,4:5)), c(1)])
  total_R_explained <- format(round(sum(varests[1:2])/sum(varests)*100,1))
  random <- rbind(random, c(total_R_explained))
  row.names(random) <- c("Intercept (ELISA_plate)",
                         "Intercept (person-level)",
                         "Residual",
                         "Number of Observations",
                         "Total R explained")
  return(random)
}

#extract fixed and random part together 
model_res.wRandomPeak <- function(model_summary){
  tTable <- res_Bchange(as.data.frame(model_summary[["tTable"]]))
  random_part <- res_randomPeak(model_summary)
  random_part[,2] <- NULL
  random_part$empty1 <- NA
  random_part$empty2 <- NA
  random_part$empty3 <- NA
  random_part$empty4 <- NA
  tTable <- rbind(tTable,c("Variance Component", rep(NA,4)))
  names(random_part) <- names(tTable)
  Res <- rbind(tTable, random_part)
  Res <- setDT(Res, keep.rownames = "predictor")
  Res$predictor[rownames(Res) == Res$predictor] <- "Random Component"
  return(Res)
}
