#------------------------------------------------------------
# 2022-05-13
# 
# Aim:
#  calculate PGSs using the metascore file
#------------------------------------------------------------
# main dir for this analysis
cd /home/ivm/analysis/PGSCatalog_Benchmarking_Sam/2022_May/Polygenic_scores


## imputed pgen files
#pgenprefix="/genesandhealth/library-red/genesandhealth/GSAv3EAMD/\
#Jul2021_44k_TOPMED-r2_Imputation_b38/topmed-r2_merged_version03/\
#chrALLincX.dose.merged_INFO0.3_MAF0.00001_F_MISSING0.2.8bit.sorted"

# I had issues when loading the pgen file directly from library-red,
# so I copied the data to a local directory as suggested by HGI.
pgenprefix="/home/ivm/analysis/GNH44k_genotype/chrALLincX.dose.merged_INFO0.3_MAF0.00001_F_MISSING0.2.8bit.sorted"
# Note - remember to DELETE THEM because storing huge data is expensive

## folder for the metascores
dir_meta="/home/ivm/analysis/PGSCatalog_Benchmarking_Sam/2022_May/\
Harmonised_metascores"




#--- plink2 using the meta-score file ---

### use a Hi_Compute VM
cd /home/ivm/analysis/PGSCatalog_Benchmarking_Sam/2022_May/Polygenic_scores/

## CVD - 1
plink2 --pfile $pgenprefix --threads 12 --memory 10000 \
--score $dir_meta/CVD_NoAmb_FirstPass.txt.gz \
header-read list-variants cols=fid,scoresums --score-col-nums 3-103 \
--out CVD_FirstPass 

## 17:33 end 18:55
## Memory usage is not high, usually <7G
## calculating AF, usually 3-6 threads; 0.5 hours 
## loading variants (can use up to 12 threads) from the meta-score file: ~ 40 min 
## finish calculating scores within a few minutes


## CVD - 2
plink2 --pfile $pgenprefix --threads 15 --memory 10000 \
--score $dir_meta/CVD_NoAmb_OtherAlleles.txt.gz \
header-read list-variants cols=fid,scoresums --score-col-nums 3-103 \
--out CVD_OtherAlleles
## 08:54 - 10:12


## Metabolic Disease - 1
plink2 --pfile $pgenprefix --threads 15 --memory 10000 \
--score $dir_meta/MetabolicDisease_NoAmb_FirstPass.txt.gz \
header-read list-variants cols=fid,scoresums --score-col-nums 3-76 \
--out MetabolicDisease_FirstPass
## 18:58 - 20:26


## Metabolic Disease - 2
plink2 --pfile $pgenprefix --threads 15 --memory 10000 \
--score $dir_meta/MetabolicDisease_NoAmb_OtherAlleles.txt.gz \
header-read list-variants cols=fid,scoresums --score-col-nums 3-76 \
--out MetabolicDisease_OtherAlleles
## 20:27 - 21:48


## Cancer - 1
plink2 --pfile $pgenprefix --threads 15 --memory 10000 \
--score $dir_meta/Cancer_NoAmb_FirstPass.txt.gz \
header-read list-variants cols=fid,scoresums --score-col-nums 3-308 \
--out Cancer_FirstPass
## 10:13 - 12:01


###### usings local pgen files
## Cancer - 2
plink2 --pfile $pgenprefix --threads 15 --memory 10000 \
--score $dir_meta/Cancer_NoAmb_OtherAlleles.txt.gz \
header-read list-variants cols=fid,scoresums --score-col-nums 3-308 \
--out Cancer_OtherAlleles
## 08:07 - 10:54


## Lipids - 1
plink2 --pfile $pgenprefix --threads 15 --memory 10000 \
--score $dir_meta/Lipids_NoAmb_FirstPass.txt.gz \
header-read list-variants cols=fid,scoresums --score-col-nums 3-54 \
--out Lipids_FirstPass
## 10:56 - 12:08 (AF calculated) - 13:23


## Lipids - 2
plink2 --pfile $pgenprefix --threads 15 --memory 10000 \
--score $dir_meta/Lipids_NoAmb_OtherAlleles.txt.gz \
header-read list-variants cols=fid,scoresums --score-col-nums 3-54 \
--out Lipids_OtherAlleles
## 14:16 - 16:43


##### using a VM with 16Gb mem and 4 CPUs
## Others - 1
plink2 --pfile $pgenprefix --threads 4 --memory 10000 \
--score $dir_meta/Others_NoAmb_FirstPass.txt.gz \
header-read list-variants cols=fid,scoresums --score-col-nums 3-264 \
--out Others_FirstPass
## 18:15 - 22:49


## Others - 2
plink2 --pfile $pgenprefix --threads 4 --memory 10000 \
--score $dir_meta/Others_NoAmb_OtherAlleles.txt.gz \
header-read list-variants cols=fid,scoresums --score-col-nums 3-264 \
--out Others_OtherAlleles
## 07:51 - 10:52










