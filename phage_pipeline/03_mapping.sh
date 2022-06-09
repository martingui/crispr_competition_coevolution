#!/bin/bash

## Arguments in the order
##      1- working directory
##      2- name of the step
##      3- name of the previous step
##      4- link to fasta reference



## DEFINE PATH
path=$1
# path=/home/enrique/work/Gandon/coevolution/phages/
step_name=$2
prev_step=$3

# ref=data/refs/Sv2972/
ref=$4



## CREATE INDEXES AND 
# bowtie2-build ${path}data/refs/Sv2972/NC_007019.1.fasta ${path}data/refs/indexes_Sv/Sv

# bowtie2-build ${path}data/refs/St/GCF_000253395.1_ASM25339v1_genomic.fna ${path}data/refs/indexes_St/St

# bowtie2 --phred33 -5 12 -p 35 -t -x ${path}data/refs/indexes_Sv/Sv -1 ${path}data/trimmed/W_seq/W4T3_S54_R1.fq.gz -2 ${path}data/trimmed/W_seq/W4T3_S54_R2.fq.gz -S ${path}results/test.sam


## CHECK FOR INDEXES

test -f  ${ref/.fasta}.1.bt2
if [ $? -eq 1 ];
    then
        echo 'Index for bowtie2 not present near fasta file'
        bowtie2-build  $ref  ${ref/.fasta}
    echo 'Created index for bowtie2'
fi





path_fasta=${path}steps/${prev_step}

path_results=${path}steps/${step_name}/

# bacteria_index=${path}data/refs/indexes_St/St

virus_index=${ref/.fasta}


for i in  $(find  $path_fasta  -name *_R1.fq.gz)
do
    ## Declare local variables
    # echo $i
    root_name=$(basename -s  _R1.fq.gz $i)
    var=$(dirname $i)
    outdir=${var/data\/trimmed/results/mapping}/

    ## Give some feedback to the user
    echo -e  "\n"phage $root_name -\> ${outdir}${root_name}.sam
    echo $i  ${i/_R1/_R2}
    echo $virus_index

    ## Mapping and indexing bam file
    echo "#### MAPPING"
    bowtie2 --phred33  -5 12  -p 24  -t  -x  $virus_index  -1 $i  -2 ${i/_R1/_R2}  -S ${outdir}${root_name}.sam 

    echo "#### SORTING"
    samtools sort  -O BAM  -o ${outdir}${root_name}.sort.bam  ${outdir}${root_name}.sam
    samtools index  -b ${outdir}${root_name}.sort.bam

done
