#! /bin/bash

## VARIABLES

# path=/home/enrique/work/Gandon/coevolution/phages/
path=$1
step_name=$2
n_threads=35
# n_threads=$3

## Maybe declaring the place where trimmomatic.jar is stored should vo to a variable



## SCRIPT

mkdir -p ${path}steps/$step_name/{W_seq,R_seq,Other_seq,uniq}

trimm_summary=${path}steps/$step_name/summary_trimm

rm -rf ${trimm_summary}_W ${trimm_summary}_R ${trimm_summary}_Other
touch ${trimm_summary}_W ${trimm_summary}_R ${trimm_summary}_Other





    ### CREATE SYMBOLIC LINKS TO SEQUENCE FILES
    ### DIIVIDED IN 3 DIRECTORIES


# echo "CREATE SYMBOLIC LINKS"

# for i in $(ls ${path}raw_data/sequences/W*); 
# do
#     ln -s $i ${path}/data/fastq_ln/W_seq/;
# done


# for i in $(ls ${path}raw_data/sequences/R*); 
# do
#     ln -s $i ${path}/data/fastq_ln/R_seq/;
# done


# for i in $(ls ${path}raw_data/sequences/ | grep -v ^W | grep -v ^R);
# do
#     ln -s ${path}raw_data/sequences/$i ${path}/data/fastq_ln/Other_seq;
# done


## DO TRIMMING FOR W_SEQUENCES.



for i in ${path}/data/fastq_ln/W_seq/*_R1_*.fastq.gz
do
    shortname=$(basename -s _L001_R1_001.fastq.gz $i);
    echo $shortname >> ${trimm_summary}_W;
    echo "Working on file " $shortname;
    java -jar /usr/local/src/Trimmomatic-0.38/trimmomatic-0.38.jar \
        PE \
        -threads $n_threads \
        -phred33 \
        -summary /tmp/tmp.trimm_summary \
        -quiet \
        $i \
        ${i/_R1_/_R2_} \
        ${path}steps/$step_name/${shortname}_R1.fq.gz \
        ${path}steps/$step_name/uniq/U1.fq.gz \
        ${path}steps/$step_name/${shortname}_R2.fq.gz \
        ${path}steps/$step_name/uniq/U2.fq.gz \
        CROP:145 \
        HEADCROP:11 \
        SLIDINGWINDOW:2:20 \
        MINLEN:50 
    cat /tmp/tmp.trimm_summary >> ${trimm_summary}_W;
    rm -rf ${path}steps/$step_name/uniq/*;
    echo -e '\n####' >> ${trimm_summary}_W;
done


for i in ${path}/data/fastq_ln/R_seq/*_R1_*.fastq.gz
do
    shortname=$(basename -s _L001_R1_001.fastq.gz $i);
    echo $shortname >> ${trimm_summary}_R;
    echo "Working on file " $shortname;
    java -jar /usr/local/src/Trimmomatic-0.38/trimmomatic-0.38.jar \
        PE \
        -threads $n_threads \
        -phred33 \
        -summary /tmp/tmp.trimm_summary \
        -quiet \
        $i \
        ${i/_R1_/_R2_} \
        ${path}steps/$step_name/R_seq/${shortname}_R1.fq.gz \
        ${path}steps/$step_name/uniq/U1.fq.gz \
        ${path}steps/$step_name/R_seq/${shortname}_R2.fq.gz \
        ${path}steps/$step_name/uniq/U2.fq.gz \
        CROP:145 \
        HEADCROP:11 \
        SLIDINGWINDOW:2:20 \
        MINLEN:50 
    cat /tmp/tmp.trimm_summary >> ${trimm_summary}_R;
    rm -rf ${path}steps/$step_name/uniq/*;
    echo -e '\n####' >> ${trimm_summary}_R;
done



for i in ${path}data/fastq_ln/Other_seq/*_R1_*.fastq.gz
do
    shortname=$(basename -s _L001_R1_001.fastq.gz $i);
    echo $shortname >> ${trimm_summary}_Other;
    echo "Working on file " $shortname;
    java -jar /usr/local/src/Trimmomatic-0.38/trimmomatic-0.38.jar \
        PE \
        -threads $n_threads \
        -phred33 \
        -summary /tmp/tmp.trimm_summary \
        -quiet \
        $i \
        ${i/_R1_/_R2_} \
        ${path}steps/$step_name/Other_seq/${shortname}_R1.fq.gz \
        ${path}steps/$step_name/uniq/U1.fq.gz \
        ${path}steps/$step_name/Other_seq/${shortname}_R2.fq.gz \
        ${path}steps/$step_name/uniq/U2.fq.gz \
        CROP:145 \
        HEADCROP:11 \
        SLIDINGWINDOW:2:20 \
        MINLEN:50 
    cat /tmp/tmp.trimm_summary >> ${trimm_summary}_Other;
    rm -rf ${path}steps/$step_name/uniq/*;
    echo -e '\n####' >> ${trimm_summary}_Other;
done

rm -rf ${path}steps/$step_name/uniq/
