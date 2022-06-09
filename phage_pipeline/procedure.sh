#! /bin/bash


## SET THE PATH TO THE PHAGES DIRECTORY
# cd /home/user/work/coev/phages 


## PREPARE DATA SUBGROUPS

## Uncompress raw data into data folder
tar -xzvf data/sequences.tar.gz -C data

## Change premissions to avoid accidents
chmod -w  data/sequences/*
chmod -w data/sequences.tar.gz

## Make links to data to make sub-groups
## It makes it easier to handle groups of files

mkdir -p  data/fastq_ln/{R_seq,W_seq,Other_seq}

for i in $PWD/data/sequences/R*
do
    ln -s  $i  $PWD/data/fastq_ln/R_seq
done


for i in $PWD/data/sequences/W*
do
    ln -s  $i  $PWD/data/fastq_ln/W_seq
done


for i in $(ls raw_data/sequences/  | grep -v -E "^W.*.fastq.gz|^R.*.fastq.gz")
do
    ln -s  $PWD/raw_data/sequences/$i  $PWD/data/fastq_ln/Other_seq
done



###################################
## CREATE PYTHON ENVIRONMENT  -- With virtualenv

./scripts/00_create_py_env.sh
## environment located in:
## ~/envs/coev/


###################################
## QUALITY CHECK 
./scripts/01_quality_check.sh  $PWD/  qc_raw_data


###################################
## Trimm and clean
./scripts/02_trimm_and_clean.sh  $PWD/  trimming


###################################
## Mapping on virus
## Arguments in the order
##      1- working directory
##      2- name of the step
##      3- name of the previous step
##      4- link to fasta reference
./scripts/03_mapping.sh  $PWD/  map_on_phage  trimming  data/refs/Sv2972/NC_007019.1_T0.fasta


###################################
## SNPcalling
## Arguments in the order
##      1- working directory
##      2- name of the step
##      3- name of the previous step. The one which produced the bam files
##      4- link to fasta reference
./scripts/04_snpcalling.sh  $PWD/  snpcall_freebayes  mapping  data/refs/Sv2972/NC_007019.1_T0.fasta 
