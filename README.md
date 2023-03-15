# Benchmarking_PGSCatalog_GnH
Calculating scores from the PGS Catalog in British Pakistani and Bangladeshi individuals in Genes &amp; Health, and assessing the performance of these scores.


## 1_scoring
We firstly calculated scores using plink2. The scoring files have been harmonised with G&H imputed variants. We received 2 multi-weight scoring files for each set of traits. Effect alleles are different for variants that are in both scoring files (to avoid flipping the effect size from the original submitted data, I think). Harmonised scoring files were prepared by Samuel Lambert. 

We next summed the two scores in R.

## 2_preparing_phenotype_table
Samples that were not inferred as British Pakistani or Bangladeshi ancestry using genetic data were excluded. We kept the unrelated subset where no pairs of individuals were inferred as 3rd degree relatives or closer. Age and sex were added to the covairate table.

## 3_validate_PGSs
R code to assess the performance of PGS Catalog scores. Curated binary phenotypes were from v005. I cleaned some quantitative phenotypes for height, BMI, WHRadjBMI, HDL, LDL, TC, TG, SBP, and DBP. Note that quatitative data have been updated since. 

- For quantitaive traits, both raw data and inverse normal transformed data were evaluated. Liner regression models were applied. Incremental R squared was reported.

- For binary traits, logistic regression models were run. Incremental AUC was reported.

In total, 791 scores were calculated and 676 of them were evaluated. Some were evaluated in a sex-specific model.
