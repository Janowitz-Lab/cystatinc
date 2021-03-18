#!/bin/bash
#$ -cwd
#$ -pe threads 12
#$ -l m_mem_free=3G
#$ -N apply_score
#$ -o score.txt
#$ -e score.txt

BINARY_UKB=/mnt/grid/ukbiobank/data/Application58510/skleeman/ukb_INFO0.8_MAF0.01_ALL
BINARY_PANIO=/mnt/grid/janowitz/rdata_norepl/pan_immuno/hail/imputed_exome_gwas_panIO_grch37

SCORE_FULL=~/PGS/final/UKB380_PGS_LDPRED2_full.tsv
SCORE_EXOME=~/PGS/final/UKB380_PGS_LDPRED2_exome.tsv

plink2 --bfile $BINARY_UKB --score $SCORE_EXOME 1 2 3 header ignore-dup-ids --allow-extra-chr --out UKB380_PGS_LDPRED2_UKB_exome


plink2 --bfile $BINARY_PANIO --score $SCORE_EXOME 4 2 3 header ignore-dup-ids --allow-extra-chr --out UKB380_PGS_LDPRED2_panIO_exome
