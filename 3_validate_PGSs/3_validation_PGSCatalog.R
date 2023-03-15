#-------------------------------------------------------
# 2023-02-13 last modified
# 
# Aims:
#  Validate scores from the PGS Catalog 
#  Quantitative data are from version 004
#  Diseases are from version 005
#-------------------------------------------------------
library(data.table)
library(foreach)
library(pROC)
library(RNOmni)
source("/home/ivm/analysis/scripts/Function_case_control_status.R")

setwd("/home/ivm/analysis/PGSCatalog_Benchmarking_Sam/2022_May/asso")


#----- Load data -----
## PGS to trait mapping file
# You will need to prepare your own file. This file has three columns: 
# (1) pgs_id: e.g. PGS002277
# (2) GnH: curated phenotypes, e.g. GNH0242_Type_2_Diabetes_narrow
# (3) Sex: "Female", "Male", and "All"; majority of the scores are tested in both sexes using all samples; some phenotypes such as breast cancer are tested in one sex; there are also some PGS that were developed for a specific sex (e.g PGS000829 and PGS000830 are two BMI PGSs developed for males and females, respectively)
PGS2trait = fread("PGS2Trait_20220620.txt")

# quantitative traits
PGS2trait_quan = PGS2trait[GnH %in% c("Height", "BMI", "HDL", "LDL", "TC", "TG",
                                      "SBP", "DBP", "WHRadjBMI")]
# binary traits
PGS2trait_binary = PGS2trait[!GnH %in% PGS2trait_quan$GnH]



## Covariate file for the unrelated subset
covunrel = fread("/home/ivm/analysis/GNH44k_phenotype/2022_06_23_Covariate_BPB_PCs_sex_age_unrelated_3rd_propIBD.txt")

## double check if all the G&H phenotype labels are correct
# GHpheno = fread("/home/ivm/analysis/GNH44k_phenotype/2022_06_21_GnH_curated_phenolist_v005.txt")
# table(unique(PGS2trait$GnH) %in% GHpheno$phenotype) 
# table(unique(PGS2trait_binary$GnH) %in% GHpheno$phenotype)  # all 50 binary traits


## quantitative data
bmi = fread("/home/ivm/analysis/GNH44k_phenotype/2022_05_26_height_weight_BMI.txt")
lipids = fread("/home/ivm/analysis/GNH44k_phenotype/2022_05_26_lipids.txt")


## PGSs
scores = fread("/home/ivm/analysis/PGSCatalog_Benchmarking_Sam/2022_May/Polygenic_scores/2022_06_29_ALL_scoresum.txt")
# keep the unrelated subset
scores = scores[OrageneID %in% covunrel$OrageneID]

#table(PGS2trait$pgs_id %in% names(scores))



#----- models -----
nPCs = 10
fullmodel = formula(paste0("pheno ~ PRS_scaled + sex + age +", 
                           paste0("PC",1:nPCs,collapse = "+")))
nullmodel = formula(paste0("pheno ~ sex + age +", 
                           paste0("PC",1:nPCs,collapse = "+")))

# sex-specific analysis
fullmodel_nosex = formula(paste0("pheno ~ PRS_scaled + age +", 
                                 paste0("PC",1:nPCs,collapse = "+")))
nullmodel_nosex = formula(paste0("pheno ~ age +", 
                                 paste0("PC",1:nPCs,collapse = "+")))




#----- test for association -----
# (1) quantitative traits
# add all quantitaive traits into one table
quancov = merge(covunrel, bmi, by = "pseudoNHSnumber", all.x = T)
quancov = merge(quancov, lipids, by = "pseudoNHSnumber", all.x = T)

# List of quantitative traits in G&H
list_quan = data.table(GnH = setdiff(names(quancov), names(covunrel)))
# PGS ID mapping to sepific GnH trait names
list_quan$label = sapply(list_quan$GnH, function(x) {
  unlist(strsplit(x, split = "_"))[1]
})
PGS2trait_quan_expand = merge(PGS2trait_quan[, .(pgs_id, label = GnH, Sex)],
                              list_quan, by = "label", allow.cartesian = T)



##----- (1) quantitative traits -----
results_quan = foreach(ii = 1:nrow(PGS2trait_quan_expand), .combine = rbind) %do% {
  myPGSID = PGS2trait_quan_expand$pgs_id[ii]
  mytrait = PGS2trait_quan_expand$GnH[ii]
  mysubgroup = PGS2trait_quan_expand$Sex[ii]
  
  cat("----- working on", mytrait, myPGSID, ";", ii, "/", nrow(PGS2trait_quan_expand) ,"-----\n")
  
  # PGS
  myPGS = scores[, .(OrageneID)]
  myPGS$PRS_scaled = scale(scores[, get(myPGSID)])
  
  # Cov and phenotype
  mydd = merge(quancov, myPGS, by = "OrageneID")
  
  # pheno
  mydd$pheno = mydd[, get(mytrait)]
  # remove missing data
  mydd = mydd[!is.na(pheno)]
  
  # sex-specific?
  if(mysubgroup == "Female") {
    mydd = mydd[sex == 2]
  } else if(mysubgroup == "Male") {
    mydd = mydd[sex == 1]
  }
  
    
  ## linear regression models
  if(mysubgroup == "All") {
    myfullmodel = fullmodel
    mynullmodel = nullmodel
  } else {
    myfullmodel = fullmodel_nosex
    mynullmodel = nullmodel_nosex
  }
  
  
  ### (1) raw data 
  myfit_full = lm(myfullmodel, data = mydd)
  myfit_CI = confint(myfit_full)
  
  myfit_null = lm(mynullmodel, data = mydd)
  
  returndd = data.table(
    GnH = mytrait, quandata = "raw", pgs_id = myPGSID, Sex = mysubgroup,
    N = nrow(mydd), prevalence = NA, Ncases = NA,
    PRS_beta = summary(myfit_full)$coefficients["PRS_scaled",1],
    PRS_beta_CIL = myfit_CI["PRS_scaled", 1], PRS_beta_CIU = myfit_CI["PRS_scaled", 2],
    PRS_se = summary(myfit_full)$coefficients["PRS_scaled",2],
    PRS_p = summary(myfit_full)$coefficients["PRS_scaled",4],
    full_R2 = summary(myfit_full)$r.squared, 
    null_R2 = summary(myfit_null)$r.squared
  )
  returndd[, incR2 := full_R2 - null_R2]
  returndd[, `:=` (AUC_full = NA, AUC_full_CIL = NA, AUC_full_CIU = NA,
                   AUC_null = NA, AUC_null_CIL = NA, AUC_null_CIU = NA, incAUC = NA)]
  
  # bootstrap to get 95% CI for R2 - takes to long
  #fullboo <- boot(mydd, function(data,indices)
  #  summary(lm(myfullmodel,data[indices,]))$r.squared,R=10000)
  #quantile(fullboo$t, c(0.025, 0.975))
  
  
  ### (2) inverse normal transformation
  mydd$pheno = RankNorm(mydd$pheno)
  
  myfit_full = lm(myfullmodel, data = mydd)
  myfit_CI = confint(myfit_full)
  
  myfit_null = lm(mynullmodel, data = mydd)
  
  returndd_INT = data.table(
    GnH = mytrait, quandata = "INT", pgs_id = myPGSID, Sex = mysubgroup,
    N = nrow(mydd), prevalence = NA, Ncases = NA,
    PRS_beta = summary(myfit_full)$coefficients["PRS_scaled",1],
    PRS_beta_CIL = myfit_CI["PRS_scaled", 1], PRS_beta_CIU = myfit_CI["PRS_scaled", 2],
    PRS_se = summary(myfit_full)$coefficients["PRS_scaled",2],
    PRS_p = summary(myfit_full)$coefficients["PRS_scaled",4],
    full_R2 = summary(myfit_full)$r.squared, 
    null_R2 = summary(myfit_null)$r.squared
  )
  returndd_INT[, incR2 := full_R2 - null_R2]
  returndd_INT[, `:=` (AUC_full = NA, AUC_full_CIL = NA, AUC_full_CIU = NA,
                   AUC_null = NA, AUC_null_CIL = NA, AUC_null_CIU = NA, incAUC = NA)]
  
  return(rbind(returndd, returndd_INT))
}


fwrite(results_quan, "2022_07_08_summary_stats_quantitative_CIs.txt", sep = "\t")



##----- (2) binary traits - no exclusion in controls -----
# diabetes - add "any diabetes"
PGS2trait_binary_diab = PGS2trait_binary[grep("Diabetes", GnH)]
anydiab = data.table(
  pgs_id = unique(PGS2trait_binary_diab$pgs_id),
  GnH = "any_Diabetes", Sex = "All"
)
PGS2trait_binary = rbind(PGS2trait_binary, anydiab)
PGS2trait_binary = PGS2trait_binary[order(pgs_id),]


results_binary = foreach(ii = 1:nrow(PGS2trait_binary), .combine = rbind) %do% {
  myPGSID = PGS2trait_binary$pgs_id[ii]
  mytrait = PGS2trait_binary$GnH[ii]
  mysubgroup = PGS2trait_binary$Sex[ii]
  # flag for diabetes
  Diabetesflag = ifelse(mytrait %in% c("any_Diabetes", "GNH0241_Type_1_Diabetes_narrow", "GNH0242_Type_2_Diabetes_narrow"), yes = T, no = F)
  
  
  cat("----- working on", mytrait, myPGSID, ";", ii, "/", nrow(PGS2trait_binary) ,"-----\n")
  
  # PGS
  myPGS = scores[, .(OrageneID)]
  myPGS$PRS_scaled = scale(scores[, get(myPGSID)])
  
  # Cov
  mydd = merge(covunrel, myPGS, by = "OrageneID")
  
  # get the phenotype
  mydd = getpheno(pheno_shortname = mytrait, 
                  pheno_type = ifelse(Diabetesflag, yes = "control exclusion", no = "binary"),
                  samplelist = mydd)
  mydd$pheno = mydd[, get(mytrait)]
  
  # remove missing data
  mydd = mydd[!is.na(pheno)]
  
  # sex-specific?
  if(mysubgroup == "Female") {
    mydd = mydd[sex == 2]
  } else if(mysubgroup == "Male") {
    mydd = mydd[sex == 1]
  }
  
  
  ## logistic regression
  if(mysubgroup == "All") {
    myfullmodel = fullmodel
    mynullmodel = nullmodel
  } else {
    myfullmodel = fullmodel_nosex
    mynullmodel = nullmodel_nosex
  }
  
  # Logistic regression
  myfit_full = glm(myfullmodel, data = mydd, family = binomial())
  myfit_CI = confint(myfit_full)
  roc_full = roc(mydd$pheno, predict(myfit_full, type = c("response")))
  
  myfit_null = glm(mynullmodel, data = mydd, family = binomial())
  roc_null = roc(mydd$pheno, predict(myfit_null, type = c("response")))
  
  returndd = data.table(
    GnH = mytrait, quandata = NA, pgs_id = myPGSID, Sex = mysubgroup,
    N = nrow(mydd), prevalence = mean(mydd$pheno), Ncases = sum(mydd$pheno),
    PRS_beta = summary(myfit_full)$coefficients["PRS_scaled",1],
    PRS_beta_CIL = myfit_CI["PRS_scaled",1],
    PRS_beta_CIU = myfit_CI["PRS_scaled",2],
    PRS_se = summary(myfit_full)$coefficients["PRS_scaled",2],
    PRS_p = summary(myfit_full)$coefficients["PRS_scaled",4],
    full_R2 = NA, null_R2 = NA, incR2 = NA,
    AUC_full = roc_full$auc[1], 
    AUC_full_CIL = ci.auc(roc_full)[1], 
    AUC_full_CIU = ci.auc(roc_full)[3],
    AUC_null = roc_null$auc[1],
    AUC_null_CIL = ci.auc(roc_null)[1], 
    AUC_null_CIU = ci.auc(roc_null)[3]
  )
  returndd[, incAUC := AUC_full - AUC_null]
  
  return(returndd)
}

# fwrite(results_binary, "2022_07_08_summary_stats_binary_CIs.txt", sep = "\t")
fwrite(results_binary, "2023_02_13_summary_stats_binary_CIs.txt", sep = "\t")



write.csv(rbind(results_quan, results_binary), "2022_07_08_summary_stats_PGSCatalog.csv", row.names = F)

