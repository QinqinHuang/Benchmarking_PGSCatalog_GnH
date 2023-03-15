#------------------------------------------------------------
# 2022-05-26
# 
# Aims: 
#  Get a list of British Pakistani and Bangladeshi people
#------------------------------------------------------------
library(data.table)
library(foreach)

setwd("~/analysis/GNH44k_phenotype")


## Teng's PCs and ancestry data
PC_ancestry = fread("/genesandhealth/red/TengHeng/FARM/GH.44k.declared_and_inferred_ethnicities.50PCs.txt")
PC_ancestry$V1 = NULL
# keep inferred P and B who are consistent with self declaration
PC_ancestry = PC_ancestry[(inferred_overall == "Bangladeshi" & declared == "1bangladeshi") | 
                            (inferred_overall == "Pakistani" & declared == "2pakistani")]
cc = paste0("PC", 1:20)
PC_BPB = cbind(PC_ancestry[, .(OrageneID, ancestry = inferred_overall)], 
               PC_ancestry[, ..cc])



## ID linkage - add pseudo NHS number
#IDlink = fread("/genesandhealth/library-red/genesandhealth/2022_03_10_pseudonhs_oragene_withmissing_DEIDENTIFIED.txt")
IDlink = fread("/genesandhealth/library-red/genesandhealth/2022_05_12_pseudoNHS_oragene_withmissing_DEIDENTIFIED.txt")
table(PC_BPB$OrageneID %in% IDlink$OrageneID)
IDlink = IDlink[!PseudoNHS %in% c("MissingNHS", "StudyRemoval", "StudyWithdrawal")]
table(PC_BPB$OrageneID %in% IDlink$OrageneID)
# 1 doesn't have a pesudoNHS number
IDlink = IDlink[OrageneID %in% PC_BPB$OrageneID]
# all IDs are unique

PC_BPB = merge(IDlink[, .(OrageneID, pseudoNHSnumber = PseudoNHS)],
               PC_BPB, by = "OrageneID")



## add year of birth and sex
QS = fread("/genesandhealth/library-red/genesandhealth/phenotypes_rawdata/QMUL__Stage1Questionnaire/2022_05_05_S1QST_redacted.tab")
QS = QS[, 1:3]
names(QS) = c("OrageneID", "MMYYYY", "sex")
# all
table(PC_BPB$OrageneID %in% QS$OrageneID)
PC_BPB_QS = merge(PC_BPB, QS, by = "OrageneID")

# reformat MM-year and calculate age at 2022-02
library(lubridate)
currtime = as.Date("01-02-2022", format = "%d-%m-%Y")
PC_BPB_QS[, age := interval(as.Date(paste0("1-", MMYYYY), format = "%d-%m-%Y"), currtime) %/% months(1)]
PC_BPB_QS[, age := age/12]


fwrite(PC_BPB_QS, "/home/ivm/analysis/GNH44k_phenotype/2022_05_26_Covariate_BPB_PCs_sex_age.txt", sep = "\t")




