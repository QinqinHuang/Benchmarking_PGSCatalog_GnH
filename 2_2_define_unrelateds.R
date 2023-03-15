#------------------------------------------------------------
# 2022-06-23
# 
# Aim:
#  make a list of unrelated samples 
#
#  -- restrict samples to those with NHS numbers (those in my
# covariate table);
#  -- merge the king table with the cov table
#  
#------------------------------------------------------------
library(data.table)

# working dir
setwd("/home/ivm/analysis/GNH44k_phenotype/")


## Load data
# relatives inferred by King PropIBD
king = fread("/home/ivm/analysis/GNH44k_genotype/GH_44k_autosome_maf0.01_geno0.01_hwe1e-6_relatedness_estimation.kin0")
# ignore 4th degree relatives and consider up to 3rd degree relatives
king = king[InfType != "4th"]

# covariate table
cov = fread("/home/ivm/analysis/GNH44k_phenotype/2022_05_26_Covariate_BPB_PCs_sex_age.txt")


## function to remove relatives
remove_rel = function(related_pairs, samplelist) {
  # relevant pairs
  related_pairs = related_pairs[(ID1 %in% samplelist$OrageneID) & (ID2 %in% samplelist$OrageneID)]
  
  cat("  Number of related pairs:", nrow(related_pairs), "\n")
  cat("   Number of unique samples:", length(c(unique(related_pairs$ID2), unique(related_pairs$ID1))), "; unique samples in col1 and col2:", length(unique(related_pairs$ID1)), length(unique(related_pairs$ID2)), "\n")
  
  # For paris with the same phenotype, remove the one with the most relatives each time
  # an empty vector to keep the list of individuals to be removed
  tbrm = c()
  # aim is to make this table empty
  leftpihat = related_pairs
  
  while(nrow(leftpihat) > 0) {
    # get the most frequent one
    rmthis = names(sort(table(c(leftpihat$ID1, leftpihat$ID2)), decreasing = T))[1]
    tbrm = c(tbrm, rmthis)
    leftpihat = leftpihat[which((!ID1 %in% rmthis) & (!ID2 %in% rmthis))]
  }
  
  cat("  ** Number of indivudals to be removed:", length(unique(tbrm)), "\n")
  samplelist = samplelist[!OrageneID %in% tbrm]
  
  return(samplelist)
}


## important to make it a character
king$ID1 = as.character(king$FID1)
king$ID2 = as.character(king$FID2)
cov$OrageneID = as.character(cov$OrageneID)

unrelated = remove_rel(king, cov)
## 29029 unrelateds

fwrite(unrelated, "2022_06_23_Covariate_BPB_PCs_sex_age_unrelated_3rd_propIBD.txt", sep = "\t")









