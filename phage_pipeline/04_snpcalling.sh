#!/bin/bash

## SNPCALLING

path=$1                     ## Working directory's path
step_name=$2
prev_step=$3                      ## name of step which created the bam files to be used
ref=$4                      ## path to reference.fasta

echo ${path}steps/
mkdir -p ${path}steps/


for i in $(find ${path}steps/${prev_step}/ -type f -name *.bam)
do
    od=$(dirname $i)
    od=${od/$prev_step/$step_name}

    mkdir -p  $od

    freebayes -b $i -f $ref -C 1 -K -F 0.01  > ${od}/$(basename -s .sort.bam  $i).vcf;

    echo $(basename -s .sort.bam  $i).vcf 'has been created';
done
