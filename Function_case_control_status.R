#------------------------------------------------------------
# 2022-06-30
# 
# Aim:
#  to make a function to prepare case-control (1-0) status
#  using the curated phenotype version 005
#
# input: 
# [pheno_shortname]: e.g. GNH0001_AtrialFibrillationFlutter
# [pheno_type]: 
#    "binary" - either a case or a control
#    "control exclusion"
# T1D 0241 narrow, exclude T2D 0242 and unspecified or rare 0244 from ctrl,
# T2D 0242 narrow, exclude T1D 0241 and unspecified or rare 0244 from ctrl
# 
# [samplelist]: a list of individuals, must have a column 
# "pseudoNHSnumber"
#
# "/home/ivm/analysis/scripts/Function_case_control_status.R"
#------------------------------------------------------------
library(data.table)

# dir for the curated phenotypes
#dirpheno = "/genesandhealth/library-red/genesandhealth/phenotypes_curated/version004_2022_02"
dirpheno = "/genesandhealth/library-red/genesandhealth/phenotypes_curated/version005_2022_06"



## function - prepare phenotype data, wrapper
getpheno = function(pheno_shortname, pheno_type, samplelist) {
  
  # (1) simple binary
  if(pheno_type == "binary") {
    return(getpheno_binary(pheno_shortname = pheno_shortname,
                           samplelist = samplelist))
    
  } else if(pheno_type == "control exclusion") {
    # (2) control exclusion
    return(getpheno_ctrlexcl(pheno_shortname = pheno_shortname,
                             samplelist = samplelist))
  }
  
}


## function - binary traits, either a case or a control,
#  this also includes T1D/T2D narrow, labelled as a complex trait
getpheno_binary = function(pheno_shortname, samplelist) {
  
  casefile = paste0(dirpheno, "/custom/Cases/", 
                    pheno_shortname, ".Cases.txt")
  
  if(!file.exists(casefile)) {
    cat(" ** Case list for", pheno_shortname, "does not exist!\n")
    return(NULL)
  } else {
    caselist = fread(casefile, header = F)
    names(caselist) = "pseudoNHSnumber"
  }
  
  statusadded = samplelist
  statusadded$newcolumn = ifelse(statusadded$pseudoNHSnumber %in% caselist$pseudoNHSnumber, yes = 1, no = 0)
  names(statusadded)[ncol(statusadded)] = pheno_shortname
  return(statusadded)
}




## function - T1D,T2D, exclude any Diabetes from ctrl
# current only accept:
# "GNH0241_Type_1_Diabetes_narrow"
# "GNH0242_Type_2_Diabetes_narrow"
# and "any_Diabetes"
getpheno_ctrlexcl = function(pheno_shortname, samplelist) {
  
  # T1D
  T1Dcase = fread(paste0(dirpheno, "/custom/Cases/GNH0241_Type_1_Diabetes_narrow.Cases.txt"), header = F)
  names(T1Dcase) = "pseudoNHSnumber"
  
  # T2D 
  T2Dcase = fread(paste0(dirpheno, "/custom/Cases/GNH0242_Type_2_Diabetes_narrow.Cases.txt"), header = F)
  names(T2Dcase) = "pseudoNHSnumber"
  
  # unspecified or rare diabetes
  otherdiab = fread(paste0(dirpheno, "/custom/Cases/GNH0244_Unspecified_or_Rare_Diabetes_narrow.Cases.txt"), header = F)
  names(otherdiab) = "pseudoNHSnumber"
  
  
  ## T1D narrow
  if(pheno_shortname == "GNH0241_Type_1_Diabetes_narrow") {
    
    statusadded = samplelist
    statusadded$newcolumn = ifelse(statusadded$pseudoNHSnumber %in% T1Dcase$pseudoNHSnumber, yes = 1, no = 0)
    
    # exclude T2D, unspecified or rare diabetes from the control set
    statusadded[statusadded$pseudoNHSnumber %in% c(T2Dcase$pseudoNHSnumber, otherdiab$pseudoNHSnumber), newcolumn := NA]
    names(statusadded)[ncol(statusadded)] = pheno_shortname
    return(statusadded)
    
    
  } else if(pheno_shortname == "GNH0242_Type_2_Diabetes_narrow") {
    ## T2D narrow
    
    statusadded = samplelist
    statusadded$newcolumn = ifelse(statusadded$pseudoNHSnumber %in% T2Dcase$pseudoNHSnumber, yes = 1, no = 0)
    
    # exclude T1D, unspecified or rare diabetes from the control set
    statusadded[statusadded$pseudoNHSnumber %in% c(T1Dcase$pseudoNHSnumber, otherdiab$pseudoNHSnumber), newcolumn := NA]
    names(statusadded)[ncol(statusadded)] = pheno_shortname
    return(statusadded)
    
    
  } else if(pheno_shortname == "any_Diabetes") {
    ## Any diabetes: T1D, T2D, or unspecified or rare diabetes
    
    statusadded = samplelist
    statusadded$newcolumn = ifelse(statusadded$pseudoNHSnumber %in% c(T1Dcase$pseudoNHSnumber, T2Dcase$pseudoNHSnumber, otherdiab$pseudoNHSnumber), yes = 1, no = 0)
    names(statusadded)[ncol(statusadded)] = pheno_shortname
    return(statusadded)
    
  } else {
    cat(" ** not processed:", pheno_shortname, "\n")
    return(NULL)
  }

}



