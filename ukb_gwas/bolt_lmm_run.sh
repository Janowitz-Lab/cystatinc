#!/bin/bash

for folder in /mnt/grid/ukbiobank/data/Application58510/skleeman/gwas_cystatinc/*/; do

base=$(basename $folder)

echo $base

cat <<EOF >./scripts/$base.sh

#!/bin/bash
#$ -cwd
#$ -pe threads 16
#$ -l m_mem_free=8G
#$ -N gwas_$base
#$ -o $base.txt
#$ -e $base.txt

cd /mnt/grid/ukbiobank/data/Application58510/skleeman/gwas_cystatinc/$base

/mnt/grid/janowitz/home/applications/BOLT-LMM_v2.3.4/bolt \
    --bfile=filtered_grch37 --noBgenIDcheck \
    --phenoFile=covariates.tsv \
    --phenoCol=egfr_cystatin \
    --covarFile=covariates.tsv \
    --covarCol=centre \
    --covarCol=array \
    --covarCol=sex \
    --covarMaxLevels=30 \
    --qCovarCol=age \
    --qCovarCol=age_squared \
    --qCovarCol=PC{1:20} \
    --LDscoresFile=ldscore.gz  --LDscoresCol=L2 \
    --geneticMapFile=/mnt/grid/janowitz/home/applications/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz \
    --lmmForceNonInf \
    --numThreads=32 \
    --statsFile=cystatin_bolt_stats.gz \
    --bgenFile=/mnt/grid/ukbiobank/data/Application58510/imputed_missing/ukb_imp_chr{1:22}_v3.bgen \
    --bgenMinMAF=0.01 \
    --bgenMinINFO=0.8 \
    --sampleFile=/mnt/grid/ukbiobank/data/Application58510/rawdata/ukb58510_imp_chr1_v3_s487282.sample \
    --statsFileBgenSnps=cystatin_bolt_stats_bgen.gz \
    --verboseStats


EOF
done
