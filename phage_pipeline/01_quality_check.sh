#!/bin/bash


## Launch example from coevolution/phages/
## ./scripts/01_quality_check.sh  $PWD/  qc_rawd_data


## PATH TO WORKING DIRECTORY
path=$1
#path=/home/enrique/work/Gandon/coevolution/phages/
step_name=$2



## START ENVIRONMENT TO EXECUTE MULTIQC

source ~/envs/coev/bin/activate

## Some temporary files are created nearby the source files
## Using links from a tmp will avoid errors in case the directory
## containing the raw data is not writable

## SEPARATE THE FILES BY NAME 

echo "CREATE TEMPORARY DIRECTORIES"

W_seq=$(mktemp -d)
R_seq=$(mktemp -d)
Other_seq=$(mktemp -d)

echo $W_seq $R_seq $Other_seq


    ### CREATE SYMBOLIC LINKS TO SEQUENCE FILES
    ### DIVIDED IN 3 DIRECTORIES

echo "CREATE SYMBOLIC LINKS"

for i in $(ls ${path}data/sequences/W*);
do
  ln -s $i $W_seq;
done


for i in $(ls ${path}data/sequences/R*);
do
  ln -s $i $R_seq;
done


for i in $(ls ${path}data/sequences/ | grep -v ^W | grep -v ^R);
do
  #echo $i;
  ln -s ${path}raw_data/sequences/$i $Other_seq
done


echo "CREATE SYMBOLIC LINKS -- DONE"

### MAKE RUN FASTQC ON EACH GROUPE

mkdir -p ${path}steps/${step_name}/fastqc/{W_seq,R_seq,Other_seq}


fastqc -t 35 --noextract -o ${path}steps/$step_name/fastqc/W_seq $W_seq/*
# multiqc -f -i W_seq -o ${path}qual/multiqc/   ${path}qual/fastqc/W_seq 

fastqc -t 35 --noextract -o ${path}steps/$step_name/fastqc/R_seq $R_seq/*
# multiqc -f -i R_seq -o ${path}qual/multiqc/ ${path}qual/fastqc/R_seq/

fastqc -t 35 --noextract -o ${path}steps/$step_name/fastqc/Other_seq $Other_seq/*
# multiqc -f -i Other_seq -o ${path}qual/multiqc/ ${path}qual/fastqc/Other_seq



## MAKE MULTIQC
# multiqc -f -i W1_seq -n W1 -o ${path}qual/multiqc/   ${path}qual/fastqc/W_seq/W1*


## LOOP MULTIQC DEPENDING ON THE INPUTS
## MAKE DIRECTORIES FOR THE DIFFERENT SAMPLES:
mkdir -p ${path}steps/$step_name/multiqc/{W_seq,R_seq,Other_seq}

for i in $(seq 8)
do
    multiqc -f -i W${i}_seq -n W${i} -o ${path}steps/$step_name/multiqc/W_seq   ${path}steps/$step_name/fastqc/W_seq/W${i}*
    multiqc -f -i R${i}_seq -n R${i} -o ${path}steps/$step_name/multiqc/R_seq   ${path}steps/$step_name/fastqc/R_seq/R${i}*
done

for i in 2972 A B C D CTRL T Undetermined
do
    multiqc -f -i ${i}_seq -n ${i} -o ${path}steps/$step_name/multiqc/Other_seq  ${path}steps/$step_name/fastqc/Other_seq/${i}*
done


rm -rf $W_seq $R_seq $Other_seq

deactivate
