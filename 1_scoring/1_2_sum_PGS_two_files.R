#------------------------------------------------------------
# 2022-06-29
# 
# Scores were calculated using two meta-score file: 
# calcualte sum
#------------------------------------------------------------
library(data.table)
library(foreach)

setwd("~/analysis/PGSCatalog_Benchmarking_Sam/2022_May/Polygenic_scores/debug")

# loop over all categories
category = c("Cancer", "CVD", "Lipids", "MetabolicDisease", "Others")

for(trait in category) {
  score1 = fread(paste0(trait, "_FirstPass.sscore"))
  names(score1)[1] = "OrageneID"
  score2 = fread(paste0(trait, "_OtherAlleles.sscore"))
  names(score2)[1] = "OrageneID"
  
  if(identical(score1$IID, score2$IID) & identical(names(score1), names(score2))) {
    s1 = as.data.frame(score1)
    s2 = as.data.frame(score2)
    s1[, -1:-2] = s1[,-1:-2] + s2[, -1:-2]
  }
  
  # keep only OraeneID
  s1$OrageneID = sapply(strsplit(s1$IID, split = "_"), function(x) {x[1]})
  
  # column names are PGSID
  names(s1) = gsub("_SUM", "", names(s1))
  
  fwrite(s1[,-2], paste0(trait, "_scoresum.txt"), sep = "\t")
}


## make a table to keep all PGSs in columns
setwd("~/analysis/PGSCatalog_Benchmarking_Sam/2022_May/Polygenic_scores/")

MD = fread("MetabolicDisease_scoresum.txt")
cancer = fread("Cancer_scoresum.txt")
CVD = fread("CVD_scoresum.txt")
lipids = fread("Lipids_scoresum.txt")
other = fread("Others_scoresum.txt")

# any duplicated PGS IDs ?
allnames = c(names(MD)[-1], names(cancer)[-1], names(CVD)[-1], names(lipids)[-1], names(other)[-1])
dup = allnames[duplicated(allnames)]

#### 4 in both CVD and MetabolicDisease
intersect(dup, names(CVD))
intersect(dup, names(MD))

PGS2trait = fread("../asso/PGS2Trait_20220620.txt")
PGS2trait[pgs_id %in% dup]
## this is for a diabetic eye disease, which is not defined in G&H
cor(CVD$PGS000819, MD$PGS000819)
cor(CVD$PGS000862, MD$PGS000862)
cor(CVD$PGS001819, MD$PGS001819)
cor(CVD$PGS002027, MD$PGS002027)
## all perfect correlation, just keep the one in the MetabolicDisease table
CVD[, `:=` (PGS000819 = NULL, PGS000862 = NULL, PGS001819 = NULL, PGS002027 = NULL)]


identical(MD$OrageneID, cancer$OrageneID)
identical(MD$OrageneID, CVD$OrageneID)
identical(MD$OrageneID, lipids$OrageneID)
identical(MD$OrageneID, other$OrageneID)

scores = cbind(MD, cancer[,-1], CVD[, -1], lipids[, -1], other[,-1])
dim(scores)
length(unique(names(scores)))
# 791 unique PGS IDs


fwrite(scores, "2022_06_29_ALL_scoresum.txt", sep = "\t")

