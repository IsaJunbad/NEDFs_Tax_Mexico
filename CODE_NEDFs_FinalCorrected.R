# Main Code for "Building upon the Nonessential Energy-Dense Foods Tax in 
#           Mexico: a modeling study of increasing the tax to improve benefits."

# Authors: Isabel Junquera-Badilla
#          Ana Basto-Abreu
#          Alan Reyes‑García
#          M. Arantxa Colchero
#          Tonatiuh Barrientos-Gutierrez

# Contact information: Center for Population Health Research, National Institute 
#                      of Public Health, Avenida Universidad 655, Santa María 
#                      Ahuacatitlán, 62100 Cuernavaca, Morelos, Mexico. 


# SETUP ########################################################################

#Clean up the environment 
rm(list=ls())      

#Set working directory to the location of the file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load required libraries 
if (!require(devtools)){install.packages("devtools")}
if (!require(bw)){devtools::install_github("INSP-RH/bw", build_vignettes = TRUE)}
if (!require(survey)){install.packages("survey")}
if (!require(tidyverse)){install.packages("tidyverse")}
if (!require(xlsx)){install.packages("xlsx")}

# DATA #########################################################################

#Import population estimates since 2019 for adults age 20 or more 
Population <- read_csv("DATA_CONAPO2018_20plus.csv")

#Read main data, BMICutoff was set to 30
DATA      <- read_csv("DATA_ENSANUT2018_NEDFs.csv") 

#Add variables needed for analysis 
years <- 10
days <- (0:years)*365
BMICutoff = 30

#Filter Population to years of interest
Population <- Population[1:years,]

#Calculate population rates comparing to sample
Population$Rate <- Population$Population_20plus/sum(DATA$Weight_svy)


# TAX SCHEMES ##################################################################

Schemes      <- c("Double", "Triple", "Quadruple")

Tax_Change   <- c(0.08, 0.16, 0.24)

Compensation <- c(0.00) 

# ANALISIS VARIABLES ###########################################################

for(j in 1:length(Tax_Change)){
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Reduction_Kcal_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Final_weight_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Reduction_weight_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Final_BMI_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Reduction_BMI_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Final_obesity_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Reduction_obesity_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("PreventedObesityCases_", Schemes[j])
  
  for (i in 1:years) {
    DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
    names(DATA)[ncol(DATA)] <- paste0("BMI_", 2018 + i, "_", Schemes[j])
    
    DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
    names(DATA)[ncol(DATA)] <- paste0("Obesity_", 2018 + i, "_", Schemes[j])
    
    DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
    names(DATA)[ncol(DATA)] <- paste0("Reduction_obesity_", 2018 + i, "_", Schemes[j])
    
    DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
    names(DATA)[ncol(DATA)] <- paste0("PreventedObesityCases_", 2018 + i, "_", Schemes[j])
  }
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Reduction_relative_obesity_", Schemes[j])
  
  for (i in unique(DATA$Sex)) {
    DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
    names(DATA)[ncol(DATA)] <- paste0("Reduction_relative_obesity_", i, "_", Schemes[j])
  }
  
  for (i in unique(DATA$SES_3t)) {
    DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
    names(DATA)[ncol(DATA)] <- paste0("Reduction_relative_obesity_SES_", i, "_", Schemes[j])
  }
  
  for (i in unique(DATA$AgeGroup)) {
    DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
    names(DATA)[ncol(DATA)] <- paste0("Reduction_relative_obesity_AgeGroup_", i, "_", Schemes[j])
  }
}

# SIMULATIONS ##################################################################

for (j in 1:length(Tax_Change)) {
  
  DATA[, paste0("Reduction_Kcal_", Schemes[j])] <- DATA$NEDF_Kcal*Tax_Change[j]*(-DATA$ChangePurchases/0.08)*(1-Compensation)
  
  #We define the matrix of daily changes in energy intake caused by the tax change
  mat_reduc    <- matrix(rep(unlist(DATA[, paste0("Reduction_Kcal_", Schemes[j])]), years + 1), ncol = years+1)
  EIchange_NEDF <- energy_build(mat_reduc, days, "Stepwise_R")
  
  rm(mat_reduc) #Eliminate mat_reduc to save memory
  
  ### We run the Hall model 
  NEDF_model <- adult_weight(bw = DATA$Weight, ht = (DATA$Height)/100,
                             age = DATA$Age, sex = DATA$Sex,
                             EIchange = EIchange_NEDF, days = max(days))
  
  rm(EIchange_NEDF) #Eliminate EIchange_NEDF to save memory
  
  ### We save anthropometric data 
  
  DATA[, paste0("Final_weight_", Schemes[j])] <- NEDF_model$Body_Weight[, 365*years ] 
  
  DATA[, paste0("Reduction_weight_", Schemes[j])] <- DATA[, paste0("Final_weight_", Schemes[j])] - DATA$Weight 
  
  DATA[, paste0("Final_BMI_", Schemes[j])]  <- NEDF_model$Body_Mass_Index[, 365*years ]
  
  DATA[, paste0("Reduction_BMI_", Schemes[j])]  <- DATA[, paste0("Final_BMI_", Schemes[j])] - DATA$BMI
  
  DATA[, paste0("Final_obesity_", Schemes[j])] <- ifelse(DATA[, paste0("Final_BMI_", Schemes[j])] >= BMICutoff, 1, 0)
  
  DATA[, paste0("Reduction_obesity_", Schemes[j])] <- DATA[, paste0("Final_obesity_", Schemes[j])] - DATA$Obesity
  
  DATA[, paste0("PreventedObesityCases_", Schemes[j])] <- DATA[, paste0("Reduction_obesity_", Schemes[j])]*as.double(Population[years, 3])
  
  
  for (i in 1:years) {
    DATA[, paste0("BMI_", 2018 + i, "_", Schemes[j])] <- NEDF_model$Body_Mass_Index[, 365 * i]
    
    DATA[, paste0("Obesity_", 2018 + i, "_", Schemes[j])]  <- ifelse(DATA[, paste0("BMI_", 2018 + i, "_", Schemes[j])] >= BMICutoff, 1, 0)
    
    DATA[, paste0("Reduction_obesity_", 2018 + i, "_", Schemes[j])] <- DATA[, paste0("Obesity_", 2018 + i, "_", Schemes[j])] - DATA$Obesity
    
    DATA[, paste0("PreventedObesityCases_", 2018 + i, "_", Schemes[j])] <- DATA[, paste0("Reduction_obesity_", 2018 + i, "_", Schemes[j])]*as.double(Population[i, 3])
  }
  
  
  DATA[, paste0("Reduction_relative_obesity_", Schemes[j])] <- DATA[, paste0("Reduction_obesity_", Schemes[j])]
  
  for (i in unique(DATA$Sex)) {
    DATA[, paste0("Reduction_relative_obesity_", i, "_", Schemes[j])] <- DATA[, paste0("Reduction_obesity_", Schemes[j])]
  }
  
  for (i in unique(DATA$SES_3t)) {
    DATA[, paste0("Reduction_relative_obesity_SES_", i, "_", Schemes[j])] <- DATA[, paste0("Reduction_obesity_", Schemes[j])]
  }
  
  for (i in unique(DATA$AgeGroup)) {
    DATA[, paste0("Reduction_relative_obesity_AgeGroup_", i, "_", Schemes[j])] <- DATA[, paste0("Reduction_obesity_", Schemes[j])]
  }
  
  
  rm(NEDF_model) #We eliminate the anthropometric simulations, since we saved the ones we need 
  
  print(Schemes[j]) #We print the scheme we are in
  
}

# SURVEY DESIGN ################################################################

#Survey design 
Design_NEDF <- svydesign(id= ~ID, strata= ~Strata, weights= ~Weight_svy, PSU= ~PSU, 
                         data= DATA, nest = TRUE)
options(survey.lonely.psu = "adjust")


#We calculate the mean of the baseline obesity prevalence and use it to make the relative change in obesity prevalence, in total and per sub-group 
mean_Obesity <- (matrix(nrow = 1, ncol = (1+length(unique(DATA$Sex))+length(unique(DATA$SES_3t))+length(unique(DATA$AgeGroup)))))

mean_Obesity[1, 1] <- svymean(~Obesity, Design_NEDF, na.rm = TRUE)
mean_Obesity[1, 2:(1+length(unique(DATA$Sex)))] <- svyby(~Obesity, ~Sex, Design_NEDF, svymean)[,2]
mean_Obesity[1, (2+length(unique(DATA$Sex))):(1+length(unique(DATA$Sex))+length(unique(DATA$SES_3t)))] <- svyby(~Obesity, ~SES_3t, Design_NEDF, svymean)[,2]
mean_Obesity[1, (2+length(unique(DATA$Sex))+length(unique(DATA$SES_3t))):ncol(mean_Obesity)] <- svyby(~Obesity, ~AgeGroup, Design_NEDF, svymean)[,2]

mean_Obesity <- data.frame(mean_Obesity)
names(mean_Obesity) <- c("Total", unique(DATA$Sex), paste0("SES_", unique(DATA$SES_3t)), paste0("AgeGroup_", unique(DATA$AgeGroup)))


for (j in 1:length(Tax_Change)) {
  
  DATA[, paste0("Reduction_relative_obesity_", Schemes[j])] <- DATA[, paste0("Reduction_relative_obesity_", Schemes[j])]/as.double(mean_Obesity[1, "Total"])
  
  for (i in unique(DATA$Sex)) {
    DATA[, paste0("Reduction_relative_obesity_", i, "_", Schemes[j])] <- DATA[, paste0("Reduction_relative_obesity_", i, "_", Schemes[j])]/as.double(mean_Obesity[1, i])
  }
  
  for (i in unique(DATA$SES_3t)) {
    DATA[, paste0("Reduction_relative_obesity_SES_", i, "_", Schemes[j])] <- DATA[, paste0("Reduction_relative_obesity_SES_", i, "_", Schemes[j])]/as.double(mean_Obesity[1, paste0("SES_", i)])
  }
  
  for (i in unique(DATA$AgeGroup)) {
    DATA[, paste0("Reduction_relative_obesity_AgeGroup_", i, "_", Schemes[j])] <- DATA[, paste0("Reduction_relative_obesity_AgeGroup_", i, "_", Schemes[j])]/as.double(mean_Obesity[1, paste0("AgeGroup_", i)])
  }
  
}

#We run the survey design with the relative change in obesity prevalence
Design_NEDF <- svydesign(id= ~ID, strata= ~Strata, weights= ~Weight_svy, PSU= ~PSU, 
                         data= DATA, nest = TRUE)
options(survey.lonely.psu = "adjust")

# POPULATION RESULTS ###########################################################

#Create new matrix to save population results
SimReduc <- (matrix(nrow = length(Tax_Change), ncol = 23))

#Fill in the matrix with the population results
for (j in 1:length(Tax_Change)) {
  SimReduc[j, 1] <- Schemes[j]
  
  SimReduc[j, 2] <- paste0(format(round(Tax_Change[j]*100, digits = 1), nsmall = 1), "%")
  
  SimReduc[j, 3] <- paste0(Compensation*100, "%")
  
  ### Calories
  Reduction_Kcal <- paste0("svymean(~Reduction_Kcal_", Schemes[j], ", Design_NEDF, na.rm = TRUE)")
  SimReduc[j, 4] <- format(round(eval(parse(text = Reduction_Kcal))[1], digits = 1), nsmall = 1)
  
  CIKCAL <- paste0("confint(", Reduction_Kcal, ")")
  SimReduc[j, 5] <- paste0("(", format(round(eval(parse(text = CIKCAL))[1], digits = 1), nsmall = 1), 
                           ", ", format(round(eval(parse(text = CIKCAL))[2], digits = 1), nsmall = 1), ")")
  
  ### Weight
  Reduction_Weight <- paste0("svymean(~Reduction_weight_", Schemes[j], ", Design_NEDF, na.rm = TRUE)")
  SimReduc[j, 6] <- format(round(eval(parse(text = Reduction_Weight))[1], digits = 1), nsmall = 1)
  
  CIW <- paste0("confint(", Reduction_Weight, ")")
  SimReduc[j, 7] <- paste0("(", format(round(eval(parse(text = CIW))[1], digits = 1), nsmall = 1),
                           ", ", format(round(eval(parse(text = CIW))[2], digits = 1), nsmall = 1), ")")
  
  ### BMI
  Reduction_BMI <- paste0("svymean(~Reduction_BMI_", Schemes[j], ", Design_NEDF, na.rm = TRUE)")
  SimReduc[j, 8] <- format(round(eval(parse(text = Reduction_BMI))[1], digits = 1), nsmall = 1)
  
  CIBMI <- paste0("confint(", Reduction_BMI, ")")
  SimReduc[j, 9] <-  paste0("(", format(round(eval(parse(text = CIBMI))[1], digits = 1), nsmall = 1),
                            ", ", format(round(eval(parse(text = CIBMI))[2], digits = 1), nsmall = 1), ")")
  
  ### Obesity
  Obesity <- paste0("svymean(~Final_obesity_", Schemes[j], ", Design_NEDF, na.rm = TRUE)")
  SimReduc[j, 10] <- format(round(eval(parse(text = Obesity))[1]*100, digits = 1), nsmall = 1)
  
  CIO <- paste0("confint(", Obesity, ")")
  SimReduc[j, 11] <- paste0("(", format(round(eval(parse(text = CIO))[1]*100, digits = 1), nsmall = 1),
                            ", ", format(round(eval(parse(text = CIO))[2]*100, digits = 1), nsmall = 1), ")")
  
  ### Reduction obesity
  Reduction_Obesity <- paste0("svymean(~Reduction_obesity_", Schemes[j], ", Design_NEDF, na.rm = TRUE)")
  SimReduc[j, 12] <- format(round(eval(parse(text = Reduction_Obesity))[1]*100, digits = 1), nsmall = 1)
  
  CIOR <- paste0("confint(", Reduction_Obesity, ")")
  SimReduc[j, 13] <-  paste0("(",  format(round(eval(parse(text = CIOR))[1]*100, digits = 1), nsmall = 1),
                             ", ", format(round(eval(parse(text = CIOR))[2]*100, digits = 1), nsmall = 1), ")")
  
  #We calculate the mean and confidence interval of initial consumption 
  SimReduc[j, 14] <- format(round(svymean(~NEDF_Kcal, Design_NEDF, na.rm = TRUE), digits = 1), nsmall = 1)
  SimReduc[j, 15] <- paste0("(", format(round(confint(svymean(~NEDF_Kcal, Design_NEDF, na.rm = TRUE))[1], digits = 1), nsmall = 1),  
                            ", ", format(round(confint(svymean(~NEDF_Kcal, Design_NEDF, na.rm = TRUE))[2], digits = 1), nsmall = 1), ")")
  
  #We calculate the mean and confidence interval of initial obesity prevalence
  SimReduc[j, 16] <- format(round(svymean(~Obesity, Design_NEDF, na.rm = TRUE)*100 , digits = 1), nsmall = 1)
  SimReduc[j, 17] <- paste0("(", format(round(confint(svymean(~Obesity, Design_NEDF, na.rm = TRUE))[1]*100, digits = 1), nsmall = 1),  
                            ", ", format(round(confint(svymean(~Obesity, Design_NEDF, na.rm = TRUE))[2]*100, digits = 1), nsmall = 1), ")")
  
  ### Obesity Cases
  ObesityCases <- paste0("svytotal(~PreventedObesityCases_", Schemes[j], ", Design_NEDF, na.rm = TRUE)")
  SimReduc[j, 18] <- format(round(eval(parse(text = ObesityCases))[1]/1000000, digits = 2), nsmall = 1)
  
  CIOC <- paste0("confint(", ObesityCases, ")")
  SimReduc[j, 19] <-  paste0("(", format(round(eval(parse(text = CIOC))[1]/1000000, digits = 2), nsmall = 1),
                             ", ", format(round(eval(parse(text = CIOC))[2]/1000000, digits = 2), nsmall = 1), ")")
  
  #We calculate the total and confidence interval of the number of individuals in Millions
  SimReduc[j, 20] <- format(round(svytotal(~Count, Design_NEDF, na.rm = TRUE)/1000000, digits = 1), nsmall = 1)
  SimReduc[j, 21] <- paste0("(", format(round(confint(svytotal(~Count, Design_NEDF, na.rm = TRUE))[1]/1000000, digits = 1), nsmall = 1),  
                            ", ", format(round(confint(svytotal(~Count, Design_NEDF, na.rm = TRUE))[2]/1000000, digits = 1), nsmall = 1), ")")
  
  ### Relative Reduction obesity
  Reduction_RelObesity <- paste0("svymean(~Reduction_relative_obesity_", Schemes[j], ", Design_NEDF, na.rm = TRUE)")
  SimReduc[j, 22] <- format(round(eval(parse(text = Reduction_RelObesity))[1]*100, digits = 1), nsmall = 1)
  
  CIORr <- paste0("confint(", Reduction_RelObesity, ")")
  SimReduc[j, 23] <- paste0("(", format(round(eval(parse(text = CIORr))[1]*100, digits = 1), nsmall = 1),
                            ", ", format(round(eval(parse(text = CIORr))[2]*100, digits = 1), nsmall = 1), ")")
  
  
  print(Schemes[j])
}

#Turn the SimReduc matrix into a data-set
SimReduc <- data.frame(SimReduc)

#Add the col names 
names(SimReduc) <- c("Scenario", "Price increase", "Compensation", 
                     "Reduction of Kcal", "Confidence interval reduction Kcal",
                     "Reduction of weight", "Confidence interval reduction weight",
                     "Reduction of BMI", "Confidence interval reduction BMI",
                     "Obesity post-tax scheme", "Confidence interval obesity post-tax scheme",
                     "Reduction of obesity", "Confidence interval reduction of obesity",
                     "Consumption pre-tax", "Confidence interval consumption pre-tax",
                     "Obesity pre-tax", "Confidence interval obesity pre-tax",
                     "Total Prevented Obesity Cases", "Confidence interval Obesity Cases",
                     "N", "Confidence interval N",
                     "Relative reduction of obesity", "Confidence interval relative reduction of obesity")

# STRATIFIED RESULTS ###########################################################

#Create new matrix to save stratified results
SimReduc_Sex <- (matrix(nrow = length(Tax_Change)*length(unique(DATA$Sex)), ncol = 23))

#Fill in the matrix with the stratified results
for (j in 1:length(Tax_Change)) {
  row_index <- (j - 1) * length(unique(DATA$Sex)) + c(1:length(unique(DATA$Sex)))
  
  #We add the schemes, tax change and compensation
  SimReduc_Sex[row_index, 1] <- paste0(Schemes[j],"_",
                                       svyby(~Count, ~Sex, Design_NEDF, svytotal)[,1])
  SimReduc_Sex[row_index, 2] <- paste0(round(Tax_Change[j]*100, digits = 1), "%")
  SimReduc_Sex[row_index, 3] <- paste0(Compensation*100, "%")
  
  ### Calories
  Reduction_Kcal <- paste0("svyby(~Reduction_Kcal_", Schemes[j],  ",~Sex, Design_NEDF, svymean)")
  SimReduc_Sex[row_index, 4] <- format(round(eval(parse(text = Reduction_Kcal))[,2], digits = 1), nsmall = 1)
  
  CIKCAL <- paste0("confint(", Reduction_Kcal, ")")
  SimReduc_Sex[row_index, 5] <- paste0("(", format(round(eval(parse(text = CIKCAL))[,1], digits = 1), nsmall = 1), 
                                       ", ", format(round(eval(parse(text = CIKCAL))[,2], digits = 1), nsmall = 1), ")")
  
  ### Weight
  Reduction_Weight <- paste0("svyby(~Reduction_weight_", Schemes[j],  ",~Sex, Design_NEDF, svymean)")
  SimReduc_Sex[row_index, 6] <- format(round(eval(parse(text = Reduction_Weight))[,2], digits = 1), nsmall = 1)
  
  CIW <- paste0("confint(", Reduction_Weight, ")")
  SimReduc_Sex[row_index, 7] <- paste0("(", format(round(eval(parse(text = CIW))[,1], digits = 1), nsmall = 1),
                                       ", ", format(round(eval(parse(text = CIW))[,2], digits = 1), nsmall = 1), ")")
  
  ### BMI
  Reduction_BMI <- paste0("svyby(~Reduction_BMI_", Schemes[j],  ",~Sex, Design_NEDF, svymean)")
  SimReduc_Sex[row_index, 8] <- format(round(eval(parse(text = Reduction_BMI))[,2], digits = 1), nsmall = 1)
  
  CIBMI <- paste0("confint(", Reduction_BMI, ")")
  SimReduc_Sex[row_index, 9] <-  paste0("(", format(round(eval(parse(text = CIBMI))[,1], digits = 1), nsmall = 1),
                                        ", ", format(round(eval(parse(text = CIBMI))[,2], digits = 1), nsmall = 1), ")")
  
  ### Obesity
  Obesity <- paste0("svyby(~Final_obesity_", Schemes[j],  ",~Sex, Design_NEDF, svymean)")
  SimReduc_Sex[row_index, 10] <- format(round(eval(parse(text = Obesity))[,2]*100, digits = 1), nsmall = 1)
  
  CIO <- paste0("confint(", Obesity, ")")
  SimReduc_Sex[row_index, 11] <- paste0("(", format(round(eval(parse(text = CIO))[,1]*100, digits = 1), nsmall = 1),
                                        ", ", format(round(eval(parse(text = CIO))[,2]*100, digits = 1), nsmall = 1), ")")
  
  ### Reduction obesity
  Reduction_Obesity <- paste0("svyby(~Reduction_obesity_", Schemes[j],  ",~Sex, Design_NEDF, svymean)")
  SimReduc_Sex[row_index, 12] <- format(round(eval(parse(text = Reduction_Obesity))[,2]*100, digits = 1), nsmall = 1)
  
  CIOR <- paste0("confint(", Reduction_Obesity, ")")
  SimReduc_Sex[row_index, 13] <-  paste0("(",  format(round(eval(parse(text = CIOR))[,1]*100, digits = 1), nsmall = 1),
                                         ", ", format(round(eval(parse(text = CIOR))[,2]*100, digits = 1), nsmall = 1), ")")
  
  #We calculate the mean and confidence interval of initial consumption 
  SimReduc_Sex[row_index, 14] <- format(round(svyby(~NEDF_Kcal, ~Sex, Design_NEDF, svymean)[,2], digits = 1), nsmall = 1)
  SimReduc_Sex[row_index, 15] <- paste0("(", format(round(confint(svyby(~NEDF_Kcal, ~Sex, Design_NEDF, svymean))[,1], digits = 1), nsmall = 1), 
                                        ", ", format(round(confint(svyby(~NEDF_Kcal, ~Sex, Design_NEDF, svymean))[,2], digits = 1), nsmall = 1), ")")
  
  #We calculate the mean and confidence interval of initial obesity prevalence
  SimReduc_Sex[row_index, 16] <- format(round(svyby(~Obesity, ~Sex, Design_NEDF, svymean)[,2]*100, digits = 1), nsmall = 1)
  SimReduc_Sex[row_index, 17] <- paste0("(", format(round(confint(svyby(~Obesity, ~Sex, Design_NEDF, svymean))[,1]*100, digits = 1), nsmall = 1), 
                                        ", ", format(round(confint(svyby(~Obesity, ~Sex, Design_NEDF, svymean))[,2]*100, digits = 1), nsmall = 1), ")")
  
  ### Obesity Cases
  ObesityCases <- paste0("svyby(~PreventedObesityCases_", Schemes[j],  ",~Sex, Design_NEDF, svytotal)")
  SimReduc_Sex[row_index, 18] <- format(round(eval(parse(text = ObesityCases))[,2]/1000000, digits = 2), nsmall = 1)
  
  CIOC <- paste0("confint(", ObesityCases, ")")
  SimReduc_Sex[row_index, 19] <-  paste0("(", format(round(eval(parse(text = CIOC))[,1]/1000000, digits = 2), nsmall = 1),
                                         ", ", format(round(eval(parse(text = CIOC))[,2]/1000000, digits = 2), nsmall = 1), ")")
  
  #We calculate the total and confidence interval of the number of individuals 
  SimReduc_Sex[row_index, 20] <- format(round(svyby(~Count, ~Sex, Design_NEDF, svytotal)[,2]/1000000, digits = 1), nsmall = 1)
  SimReduc_Sex[row_index, 21] <- paste0("(", format(round(confint(svyby(~Count, ~Sex, Design_NEDF, svytotal))[,1]/1000000, digits = 1), nsmall = 1), 
                                        ", ", format(round(confint(svyby(~Count, ~Sex, Design_NEDF, svytotal))[,2]/1000000, digits = 1), nsmall = 1), ")")
  
  ### Relative Reduction obesity
  Reduction_RelObesity <- paste0("svyby(~Reduction_relative_obesity_", unique(DATA$Sex), "_", Schemes[j],  ",~Sex, Design_NEDF, svymean)")
  SimReduc_Sex[row_index, 22] <- sapply(seq_along(Reduction_RelObesity), function(i) {
                                        format(round(eval(parse(text = Reduction_RelObesity[i]))[i,2]*100, digits = 1), nsmall = 1)})
  
  CIORr <- paste0("confint(", Reduction_RelObesity, ")")
  SimReduc_Sex[row_index, 23] <-  sapply(seq_along(CIORr), function(i) {
                                         paste0("(", format(round(eval(parse(text = CIORr[i]))[i,1]*100, digits = 1), nsmall = 1),
                                               ", ", format(round(eval(parse(text = CIORr[i]))[i,2]*100, digits = 1), nsmall = 1), ")")})
  print(Schemes[j])
}

#Turn the SimReduc matrix into a data-set
SimReduc_Sex <- data.frame(SimReduc_Sex)

#Add the col names 
names(SimReduc_Sex) <- c("Scenario and Group", "Price increase", "Compensation", 
                         "Reduction of Kcal", "Confidence interval reduction Kcal",
                         "Reduction of weight", "Confidence interval reduction weight",
                         "Reduction of BMI", "Confidence interval reduction BMI",
                         "Obesity post-tax scheme", "Confidence interval obesity post-tax scheme",
                         "Reduction of obesity", "Confidence interval reduction of obesity",
                         "Consumption pre-tax", "Confidence interval consumption pre-tax",
                         "Obesity pre-tax", "Confidence interval obesity pre-tax",
                         "Total Prevented Obesity Cases", "Confidence interval Obesity Cases",
                         "N", "Confidence interval N",
                         "Relative reduction of obesity", "Confidence interval relative reduction of obesity")

#Create new matrix to save stratified results
SimReduc_SES_3t <- (matrix(nrow = length(Tax_Change)*length(unique(DATA$SES_3t)), ncol = 23))

#Fill in the matrix with the stratified results
for (j in 1:length(Tax_Change)) {
  row_index <- (j - 1) * length(unique(DATA$SES_3t)) + c(1:length(unique(DATA$SES_3t)))
  
  #We add the schemes, tax change and compensation
  SimReduc_SES_3t[row_index, 1] <- paste0(Schemes[j], "_SES ",
                                          svyby(~Count, ~SES_3t, Design_NEDF, svytotal)[,1])
  SimReduc_SES_3t[row_index, 2] <- paste0(round(Tax_Change[j]*100, digits = 1), "%")
  SimReduc_SES_3t[row_index, 3] <- paste0(Compensation*100, "%")
  
  ### Calories
  Reduction_Kcal <- paste0("svyby(~Reduction_Kcal_", Schemes[j],  ",~SES_3t, Design_NEDF, svymean)")
  SimReduc_SES_3t[row_index, 4] <- format(round(eval(parse(text = Reduction_Kcal))[,2], digits = 1), nsmall = 1)
  
  CIKCAL <- paste0("confint(", Reduction_Kcal, ")")
  SimReduc_SES_3t[row_index, 5] <- paste0("(", format(round(eval(parse(text = CIKCAL))[,1], digits = 1), nsmall = 1), 
                                          ", ", format(round(eval(parse(text = CIKCAL))[,2], digits = 1), nsmall = 1), ")")
  
  ### Weight
  Reduction_Weight <- paste0("svyby(~Reduction_weight_", Schemes[j],  ",~SES_3t, Design_NEDF, svymean)")
  SimReduc_SES_3t[row_index, 6] <- format(round(eval(parse(text = Reduction_Weight))[,2], digits = 1), nsmall = 1)
  
  CIW <- paste0("confint(", Reduction_Weight, ")")
  SimReduc_SES_3t[row_index, 7] <- paste0("(", format(round(eval(parse(text = CIW))[,1], digits = 1), nsmall = 1),
                                          ", ", format(round(eval(parse(text = CIW))[,2], digits = 1), nsmall = 1), ")")
  
  ### BMI
  Reduction_BMI <- paste0("svyby(~Reduction_BMI_", Schemes[j],  ",~SES_3t, Design_NEDF, svymean)")
  SimReduc_SES_3t[row_index, 8] <- format(round(eval(parse(text = Reduction_BMI))[,2], digits = 1), nsmall = 1)
  
  CIBMI <- paste0("confint(", Reduction_BMI, ")")
  SimReduc_SES_3t[row_index, 9] <-  paste0("(", format(round(eval(parse(text = CIBMI))[,1], digits = 1), nsmall = 1),
                                           ", ", format(round(eval(parse(text = CIBMI))[,2], digits = 1), nsmall = 1), ")")
  
  ### Obesity
  Obesity <- paste0("svyby(~Final_obesity_", Schemes[j],  ",~SES_3t, Design_NEDF, svymean)")
  SimReduc_SES_3t[row_index, 10] <- format(round(eval(parse(text = Obesity))[,2]*100, digits = 1), nsmall = 1)
  
  CIO <- paste0("confint(", Obesity, ")")
  SimReduc_SES_3t[row_index, 11] <- paste0("(", format(round(eval(parse(text = CIO))[,1]*100, digits = 1), nsmall = 1),
                                           ", ", format(round(eval(parse(text = CIO))[,2]*100, digits = 1), nsmall = 1), ")")
  
  ### Reduction obesity
  Reduction_Obesity <- paste0("svyby(~Reduction_obesity_", Schemes[j],  ",~SES_3t, Design_NEDF, svymean)")
  SimReduc_SES_3t[row_index, 12] <- format(round(eval(parse(text = Reduction_Obesity))[,2]*100, digits = 1), nsmall = 1)
  
  CIOR <- paste0("confint(", Reduction_Obesity, ")")
  SimReduc_SES_3t[row_index, 13] <-  paste0("(",  format(round(eval(parse(text = CIOR))[,1]*100, digits = 1), nsmall = 1),
                                            ", ", format(round(eval(parse(text = CIOR))[,2]*100, digits = 1), nsmall = 1), ")")
  
  #We calculate the mean and confidence interval of initial consumption 
  SimReduc_SES_3t[row_index, 14] <- format(round(svyby(~NEDF_Kcal, ~SES_3t, Design_NEDF, svymean)[,2], digits = 1), nsmall = 1)
  SimReduc_SES_3t[row_index, 15] <- paste0("(", format(round(confint(svyby(~NEDF_Kcal, ~SES_3t, Design_NEDF, svymean))[,1], digits = 1), nsmall = 1), 
                                           ", ", format(round(confint(svyby(~NEDF_Kcal, ~SES_3t, Design_NEDF, svymean))[,2], digits = 1), nsmall = 1), ")")
  
  #We calculate the mean and confidence interval of initial obesity prevalence
  SimReduc_SES_3t[row_index, 16] <- format(round(svyby(~Obesity, ~SES_3t, Design_NEDF, svymean)[,2]*100, digits = 1), nsmall = 1)
  SimReduc_SES_3t[row_index, 17] <- paste0("(", format(round(confint(svyby(~Obesity, ~SES_3t, Design_NEDF, svymean))[,1]*100, digits = 1), nsmall = 1), 
                                           ", ", format(round(confint(svyby(~Obesity, ~SES_3t, Design_NEDF, svymean))[,2]*100, digits = 1), nsmall = 1), ")")
  
  ### Obesity Cases
  ObesityCases <- paste0("svyby(~PreventedObesityCases_", Schemes[j],  ",~SES_3t, Design_NEDF, svytotal)")
  SimReduc_SES_3t[row_index, 18] <- format(round(eval(parse(text = ObesityCases))[,2]/1000000, digits = 2), nsmall = 1)
  
  CIOC <- paste0("confint(", ObesityCases, ")")
  SimReduc_SES_3t[row_index, 19] <-  paste0("(", format(round(eval(parse(text = CIOC))[,1]/1000000, digits = 2), nsmall = 1),
                                            ", ", format(round(eval(parse(text = CIOC))[,2]/1000000, digits = 2), nsmall = 1), ")")
  
  #We calculate the total and confidence interval of the number of individuals in Millions
  SimReduc_SES_3t[row_index, 20] <- format(round(svyby(~Count, ~SES_3t, Design_NEDF, svytotal)[,2]/1000000, digits = 1), nsmall = 1)
  SimReduc_SES_3t[row_index, 21] <- paste0("(", format(round(confint(svyby(~Count, ~SES_3t, Design_NEDF, svytotal))[,1]/1000000, digits = 1), nsmall = 1), 
                                           ", ", format(round(confint(svyby(~Count, ~SES_3t, Design_NEDF, svytotal))[,2]/1000000, digits = 1), nsmall = 1), ")")
  
  ### Relative Reduction obesity
  Reduction_RelObesity <- paste0("svyby(~Reduction_relative_obesity_SES_", unique(DATA$SES_3t), "_", Schemes[j],  ",~SES_3t, Design_NEDF, svymean)")
  SimReduc_SES_3t[row_index, 22] <- sapply(seq_along(Reduction_RelObesity), function(i) {
                                           format(round(eval(parse(text = Reduction_RelObesity[i]))[i,2]*100, digits = 1), nsmall = 1)})
  
  CIORr <- paste0("confint(", Reduction_RelObesity, ")")
  SimReduc_SES_3t[row_index, 23] <-  sapply(seq_along(CIORr), function(i) {
                                            paste0("(", format(round(eval(parse(text = CIORr[i]))[i,1]*100, digits = 1), nsmall = 1),
                                                  ", ", format(round(eval(parse(text = CIORr[i]))[i,2]*100, digits = 1), nsmall = 1), ")")})
                                          
  print(Schemes[j])
}

#Turn the SimReduc matrix into a data-set
SimReduc_SES_3t <- data.frame(SimReduc_SES_3t)

#Add the col names 
names(SimReduc_SES_3t) <- c("Scenario and Group", "Price increase", "Compensation",
                            "Reduction of Kcal", "Confidence interval reduction Kcal",
                            "Reduction of weight", "Confidence interval reduction weight",
                            "Reduction of BMI", "Confidence interval reduction BMI",
                            "Obesity post-tax scheme", "Confidence interval obesity post-tax scheme",
                            "Reduction of obesity", "Confidence interval reduction of obesity",
                            "Consumption pre-tax", "Confidence interval consumption pre-tax",
                            "Obesity pre-tax", "Confidence interval obesity pre-tax",
                            "Total Prevented Obesity Cases", "Confidence interval Obesity Cases",
                            "N", "Confidence interval N",
                            "Relative reduction of obesity", "Confidence interval relative reduction of obesity")

#Create new matrix to save stratified results
SimReduc_AgeGroup <- (matrix(nrow = length(Tax_Change)*length(unique(DATA$AgeGroup)), ncol = 23))

#Fill in the matrix with the stratified results
for (j in 1:length(Tax_Change)) {
  row_index <- (j - 1) * length(unique(DATA$AgeGroup)) + c(1:length(unique(DATA$AgeGroup)))
  
  #We add the schemes, tax change and compensation
  SimReduc_AgeGroup[row_index, 1] <- paste0(Schemes[j], "_AgeGroup ",
                                            svyby(~Count, ~AgeGroup, Design_NEDF, svytotal)[,1])
  SimReduc_AgeGroup[row_index, 2] <- paste0(round(Tax_Change[j]*100, digits = 1), "%")
  SimReduc_AgeGroup[row_index, 3] <- paste0(Compensation*100, "%")
  
  ### Calories
  Reduction_Kcal <- paste0("svyby(~Reduction_Kcal_", Schemes[j],  ",~AgeGroup, Design_NEDF, svymean)")
  SimReduc_AgeGroup[row_index, 4] <- format(round(eval(parse(text = Reduction_Kcal))[,2], digits = 1), nsmall = 1)
  
  CIKCAL <- paste0("confint(", Reduction_Kcal, ")")
  SimReduc_AgeGroup[row_index, 5] <- paste0("(", format(round(eval(parse(text = CIKCAL))[,1], digits = 1), nsmall = 1), 
                                            ", ", format(round(eval(parse(text = CIKCAL))[,2], digits = 1), nsmall = 1), ")")
  
  ### Weight
  Reduction_Weight <- paste0("svyby(~Reduction_weight_", Schemes[j],  ",~AgeGroup, Design_NEDF, svymean)")
  SimReduc_AgeGroup[row_index, 6] <- format(round(eval(parse(text = Reduction_Weight))[,2], digits = 1), nsmall = 1)
  
  CIW <- paste0("confint(", Reduction_Weight, ")")
  SimReduc_AgeGroup[row_index, 7] <- paste0("(", format(round(eval(parse(text = CIW))[,1], digits = 1), nsmall = 1),
                                            ", ", format(round(eval(parse(text = CIW))[,2], digits = 1), nsmall = 1), ")")
  
  ### BMI
  Reduction_BMI <- paste0("svyby(~Reduction_BMI_", Schemes[j],  ",~AgeGroup, Design_NEDF, svymean)")
  SimReduc_AgeGroup[row_index, 8] <- format(round(eval(parse(text = Reduction_BMI))[,2], digits = 1), nsmall = 1)
  
  CIBMI <- paste0("confint(", Reduction_BMI, ")")
  SimReduc_AgeGroup[row_index, 9] <-  paste0("(", format(round(eval(parse(text = CIBMI))[,1], digits = 1), nsmall = 1),
                                             ", ", format(round(eval(parse(text = CIBMI))[,2], digits = 1), nsmall = 1), ")")
  
  ### Obesity
  Obesity <- paste0("svyby(~Final_obesity_", Schemes[j],  ",~AgeGroup, Design_NEDF, svymean)")
  SimReduc_AgeGroup[row_index, 10] <- format(round(eval(parse(text = Obesity))[,2]*100, digits = 1), nsmall = 1)
  
  CIO <- paste0("confint(", Obesity, ")")
  SimReduc_AgeGroup[row_index, 11] <- paste0("(", format(round(eval(parse(text = CIO))[,1]*100, digits = 1), nsmall = 1),
                                             ", ", format(round(eval(parse(text = CIO))[,2]*100, digits = 1), nsmall = 1), ")")
  
  ### Reduction obesity
  Reduction_Obesity <- paste0("svyby(~Reduction_obesity_", Schemes[j],  ",~AgeGroup, Design_NEDF, svymean)")
  SimReduc_AgeGroup[row_index, 12] <- format(round(eval(parse(text = Reduction_Obesity))[,2]*100, digits = 1), nsmall = 1)
  
  CIOR <- paste0("confint(", Reduction_Obesity, ")")
  SimReduc_AgeGroup[row_index, 13] <-  paste0("(",  format(round(eval(parse(text = CIOR))[,1]*100, digits = 1), nsmall = 1),
                                              ", ", format(round(eval(parse(text = CIOR))[,2]*100, digits = 1), nsmall = 1), ")")
  
  #We calculate the mean and confidence interval of initial consumption 
  SimReduc_AgeGroup[row_index, 14] <- format(round(svyby(~NEDF_Kcal, ~AgeGroup, Design_NEDF, svymean)[,2], digits = 1), nsmall = 1)
  SimReduc_AgeGroup[row_index, 15] <- paste0("(", format(round(confint(svyby(~NEDF_Kcal, ~AgeGroup, Design_NEDF, svymean))[,1], digits = 1), nsmall = 1), 
                                             ", ", format(round(confint(svyby(~NEDF_Kcal, ~AgeGroup, Design_NEDF, svymean))[,2], digits = 1), nsmall = 1), ")")
  
  #We calculate the mean and confidence interval of initial obesity prevalence
  SimReduc_AgeGroup[row_index, 16] <- format(round(svyby(~Obesity, ~AgeGroup, Design_NEDF, svymean)[,2]*100, digits = 1), nsmall = 1)
  SimReduc_AgeGroup[row_index, 17] <- paste0("(", format(round(confint(svyby(~Obesity, ~AgeGroup, Design_NEDF, svymean))[,1]*100, digits = 1), nsmall = 1), 
                                             ", ", format(round(confint(svyby(~Obesity, ~AgeGroup, Design_NEDF, svymean))[,2]*100, digits = 1), nsmall = 1), ")")
  
  ### Obesity Cases
  ObesityCases <- paste0("svyby(~PreventedObesityCases_", Schemes[j],  ",~AgeGroup, Design_NEDF, svytotal)")
  SimReduc_AgeGroup[row_index, 18] <- format(round(eval(parse(text = ObesityCases))[,2]/1000000, digits = 2), nsmall = 1)
  
  CIOC <- paste0("confint(", ObesityCases, ")")
  SimReduc_AgeGroup[row_index, 19] <-  paste0("(", format(round(eval(parse(text = CIOC))[,1]/1000000, digits = 2), nsmall = 1),
                                              ", ", format(round(eval(parse(text = CIOC))[,2]/1000000, digits = 2), nsmall = 1), ")")
  
  #We calculate the total and confidence interval of the number of individuals in Millions
  SimReduc_AgeGroup[row_index, 20] <- format(round(svyby(~Count, ~AgeGroup, Design_NEDF, svytotal)[,2]/1000000, digits = 1), nsmall = 1)
  SimReduc_AgeGroup[row_index, 21] <- paste0("(", format(round(confint(svyby(~Count, ~AgeGroup, Design_NEDF, svytotal))[,1]/1000000, digits = 1), nsmall = 1), 
                                             ", ", format(round(confint(svyby(~Count, ~AgeGroup, Design_NEDF, svytotal))[,2]/1000000, digits = 1), nsmall = 1), ")")
  
  ### Relative Reduction obesity
  Reduction_RelObesity <- paste0("svyby(~Reduction_relative_obesity_AgeGroup_", unique(DATA$AgeGroup), "_", Schemes[j],  ",~AgeGroup, Design_NEDF, svymean)")
  SimReduc_AgeGroup[row_index, 22] <- sapply(seq_along(Reduction_RelObesity), function(i) {
                                             format(round(eval(parse(text = Reduction_RelObesity[i]))[i,2]*100, digits = 1), nsmall = 1)})
  
  CIORr <- paste0("confint(", Reduction_RelObesity, ")")
  SimReduc_AgeGroup[row_index, 23] <-  sapply(seq_along(CIORr), function(i) {
                                              paste0("(", format(round(eval(parse(text = CIORr[i]))[i,1]*100, digits = 1), nsmall = 1),
                                                    ", ", format(round(eval(parse(text = CIORr[i]))[i,2]*100, digits = 1), nsmall = 1), ")")})
  
  print(Schemes[j])
}

#Turn the SimReduc matrix into a data-set
SimReduc_AgeGroup <- data.frame(SimReduc_AgeGroup)

#Add the col names 
names(SimReduc_AgeGroup) <- c("Scenario and Group", "Price increase", "Compensation",
                              "Reduction of Kcal", "Confidence interval reduction Kcal",
                              "Reduction of weight", "Confidence interval reduction weight",
                              "Reduction of BMI", "Confidence interval reduction BMI",
                              "Obesity post-tax scheme", "Confidence interval obesity post-tax scheme",
                              "Reduction of obesity", "Confidence interval reduction of obesity",
                              "Consumption pre-tax", "Confidence interval consumption pre-tax",
                              "Obesity pre-tax", "Confidence interval obesity pre-tax",
                              "Total Prevented Obesity Cases", "Confidence interval Obesity Cases",
                              "N", "Confidence interval N",
                              "Relative reduction of obesity", "Confidence interval relative reduction of obesity")

# RESULTS FOR ALL YEARS ########################################################

#Create new matrix to save population results
SimReduc_Years <- (matrix(nrow = length(Tax_Change)*3, ncol = years+4))

#Fill in the matrix with the population results
for (j in 1:length(Tax_Change)) {
  
  row_index <- ((1:(length(Tax_Change)*3))[(1:(length(Tax_Change)*3)) %% 3 == 1])[j]
  
  SimReduc_Years[row_index  , 1] <- paste0(Schemes[j], "_ObesityPrevalence")
  SimReduc_Years[row_index+1, 1] <- paste0(Schemes[j], "_ChangeInPrevalence")
  SimReduc_Years[row_index+2, 1] <- paste0(Schemes[j], "_PreventedObesityCases")
  
  SimReduc_Years[row_index:(row_index+2), 2] <- paste0(format(round((Tax_Change[j]+0.08)*100, digits = 1), nsmall = 1), "%")
  SimReduc_Years[row_index:(row_index+2), 3] <- paste0(Compensation*100, "%")
  
  SimReduc_Years[row_index  , 4] <- format(round(svymean(~Obesity, Design_NEDF, na.rm = TRUE)*100 , digits = 10), nsmall = 1)
  SimReduc_Years[row_index+1, 4] <- 0
  SimReduc_Years[row_index+2, 4] <- 0
  
  for (i in 1:years) {
    Obesity <- paste0("svymean(~Obesity_", 2018 + i, "_", Schemes[j], ", Design_NEDF, na.rm = TRUE)")
    SimReduc_Years[row_index  , i+4] <- round(eval(parse(text = Obesity))[1]*100, digits = 10)
    
    Reduction_Obesity <- paste0("svymean(~Reduction_obesity_", 2018 + i, "_", Schemes[j], ", Design_NEDF, na.rm = TRUE)")
    SimReduc_Years[row_index+1, i+4] <- round(eval(parse(text = Reduction_Obesity))[1]*100, digits = 10)
    
    ObesityCases <- paste0("svytotal(~PreventedObesityCases_", 2018 + i, "_", Schemes[j], ", Design_NEDF, na.rm = TRUE)")
    SimReduc_Years[row_index+2, i+4] <- round(eval(parse(text = ObesityCases))[1], digits = 10)
  }  
  
  print(Schemes[j])
  
}

#Turn the SimReduc matrix into a data-set
SimReduc_Years <- data.frame(SimReduc_Years)

#Add the col names 
names(SimReduc_Years) <- c("Scenario and Outcome", "NEDFs Tax", "Compensation", paste("Year",2018:(2018+years)))


# SAVE ALL RESULTS #############################################################

name_exportFile = paste0("NEDF 2018_",Sys.Date(),"_","FinalCorrected.xlsx")

# Write the first data set in a new workbook
write.xlsx(SimReduc,
           row.names=FALSE,
           file = name_exportFile, 
           sheetName = "Population Results", 
           append = FALSE)

# Add the other data sets in new worksheets
write.xlsx(SimReduc_Sex,  
           row.names=FALSE,
           file = name_exportFile, 
           sheetName="Stratified by Sex", 
           append=TRUE)

write.xlsx(SimReduc_SES_3t,  
           row.names=FALSE,
           file = name_exportFile, 
           sheetName="Stratified by Terciles of income", 
           append=TRUE)

write.xlsx(SimReduc_AgeGroup,  
           row.names=FALSE,
           file = name_exportFile, 
           sheetName="Stratified by Age groups", 
           append=TRUE)

write.xlsx(SimReduc_Years,  
           row.names=FALSE,
           file = name_exportFile, 
           sheetName="Population Results by year", 
           append=TRUE)

beepr::beep(4)

# PLOTS ########################################################################

# Pivot to long format
SimReduc_Years_long <- SimReduc_Years %>%
  pivot_longer(cols = starts_with("Year"), names_to = "Year", values_to = "Value") %>%
  mutate(Year = as.numeric(sub("Year ", "", Year))) %>%
  mutate(Value = as.numeric(Value))

# Modify the Scenario and Outcome column to be a factor with specified levels
SimReduc_Years_long <- SimReduc_Years_long %>%
  mutate(`Scenario and Outcome` = factor(`Scenario and Outcome`, levels = unique(SimReduc_Years$`Scenario and Outcome`)))%>%
  mutate(Outcome = sub(".*_", "", `Scenario and Outcome`))

# Split the data into three separate data sets based on the outcome
SimReduc_ObesityPrevalence <- SimReduc_Years_long %>%
  filter(Outcome == "ObesityPrevalence")

SimReduc_ChangeInPrevalence <- SimReduc_Years_long %>%
  filter(Outcome == "ChangeInPrevalence")

SimReduc_PreventedObesityCases <- SimReduc_Years_long %>%
  filter(Outcome == "PreventedObesityCases")


#Plot population trends of the results 
plot_ObesityPrevalence <- ggplot(data = SimReduc_ObesityPrevalence, aes(x = Year, y = Value, color = `NEDFs Tax`)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Obesity Prevalence trend",
       x = "Year", y = "%",
       color = "NEDFs Tax") +
  scale_x_continuous(breaks = seq(min(SimReduc_ObesityPrevalence$Year), max(SimReduc_ObesityPrevalence$Year), by = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
  theme_bw(); plot_ObesityPrevalence

plot_ChangeInPrevalence <- ggplot(data = SimReduc_ChangeInPrevalence, aes(x = Year, y = Value, color = `NEDFs Tax`)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Change in Obesity Prevalence",
       x = "Year", y = "pp",
       color = "NEDFs Tax") +
  scale_x_continuous(breaks = seq(min(SimReduc_ChangeInPrevalence$Year), max(SimReduc_ChangeInPrevalence$Year), by = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
  theme_bw(); plot_ChangeInPrevalence

plot_PreventedObesityCases <- ggplot(data = SimReduc_PreventedObesityCases, aes(x = Year, y = Value, color = `NEDFs Tax`)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Prevented Obesity Cases",
       x = "Year", y = "Number of people",
       color = "NEDFs Tax") +
  scale_x_continuous(breaks = seq(min(SimReduc_PreventedObesityCases$Year), max(SimReduc_PreventedObesityCases$Year), by = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
  theme_bw(); plot_PreventedObesityCases

#Save plots
ggsave( filename = paste0("PLOT_ObesityPrevalenceTrend.png"),
        plot = plot_ObesityPrevalence, width = 13, height = 4.3, dpi = 300)
ggsave( filename = paste0("PLOT_ChangeObesityPrevalence.png"),
        plot = plot_ChangeInPrevalence, width = 13, height = 4.3, dpi = 300)
ggsave( filename = paste0("PLOT_PreventesCasas.png"),
        plot = plot_PreventedObesityCases, width = 13, height = 4.3, dpi = 300)
