# Polygenic risk score
- Construct_PRS.ipynb: Jupyter Notebook (R kernel) with documented code to go from GWAS summary statistics to internally tuned polygenic risk score
- score_plink.sh: PLINK2 script for linear scoring using PRS generated from UKB training set, implemented in LDPred2
- UKB380_PGS_LDPRED2.tsv: Weights for cystatin C production polygenic risk score trained and tuned in 380k european subjects from UKB, uses high-quality HapMap3 SNPs provided by authors
- UKB380_PGS_LDPRED2_exome.tsv: As above but on the subset of SNPs (around 300k) that are captured in the imputed exome sequencing data.
