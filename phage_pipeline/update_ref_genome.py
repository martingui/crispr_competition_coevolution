#!/usr/bin/python3

from Bio import SeqIO
from scripts import vcf_parser3 as vp
import numpy as np

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


## This script allows to create a consensus sequence froma
## FASTA file and a VCF file
## The equivalent in bcftools is:
## bcftools consensus --sample unknown -f NC_007019.1.fasta TO-WT_S83.vcf.gz -o test.fasta



## Input fasta file
fasta_file = "data/NC_007019.1.fasta"

## Output file
outfile = fasta_file.replace('.fasta', '_T0.fasta').replace('data','steps/new_ref_T0')

## Minimum haplotype frequency to replace REF genotype with ALT genotype
min_hap_freq = 0.45
print("\nOnly mutations with an Allele Frequency higher than {} will be retained\n".format(min_hap_freq))



## Import sequence
ref = SeqIO.read(fasta_file, "fasta")


## Copy the sequence as a string into a variable
seq_str = str(ref.seq)


## Import vcf to be used as a guide for replacements
vcf_file = vp.ImportVCF2()
vcf_file.load_vcf('data/TO-WT_S83.vcf')

counter = 0

## Loop through the reversed list of indexes.
for i in reversed(vcf_file.index): 
    ## Create a list containing the Allele Frequency for that Variant
    my_list = list(np.array(vcf_file.INFO[i]['AO'])/vcf_file.INFO[i]['DP'])

    ## Extra info for the stdout
    # print("Mutation ID= {} ; AF= {}".format(i, my_list))

    ## Filter out Variants with a low Allele Frequency.
    ## Only keep the Most frequent Variant
    if max(my_list) > min_hap_freq:

        ## Extra info for the stdout
        # print("Retained Mutation")

        counter +=1
        ## Get the index of the Variant (haplotype) 
        hap_index = my_list.index(max(my_list))
        ## Get the position|locus for the variant with the higher AF.
        ## Remove 1 to make it fit with python indexes of strings
        pos_str = vcf_file.POS[i] -1
        ## Get the lenght of the Variant
        length = len(vcf_file.REF[i])
        ## Create a new string with the sequence Before the variant + Alt variant + sequence After variant
        seq_str = seq_str[:pos_str] + str(vcf_file.ALT[i][hap_index]) + seq_str[pos_str + length :]

        ## Print stuff fot the user
        print("DP= {}, AO= {}".format(vcf_file.INFO[i]['DP'], my_list))
        print( "Index= {} ; ALT= {} ; REF= {} ; POS= {} ".format(hap_index, vcf_file.ALT[i][hap_index], vcf_file.REF[i], vcf_file.POS[i]) )
        print( "POS= {} ; REF_found={} ;".format(pos_str+1, seq_str[pos_str:pos_str+length] ))
        print("REF= " + seq_str[pos_str-5:pos_str].lower() + \
            str(vcf_file.REF[i]) + \
            seq_str[pos_str + length:pos_str + length + 4].lower() )
        print("ALT= " + seq_str[pos_str-5:pos_str].lower() + str(vcf_file.ALT[i][hap_index]) + seq_str[pos_str + length : pos_str + length + 4].lower() + '\n')

    ## Extra info for the stdout
    # else: 
        # print("Not retained\n")


print("Original length of sequence: {}".format(len(ref)))
print("New sequence length : {} \n".format(len(seq_str)))

print("Number of integrated Mutations : {}\n".format(counter))

## Convert sequence string into a sequence
seq_object = Seq(seq_str, 'generic_dna')

## Replace attributes in original sequence object
ref.seq = seq_object
ref.id = ref.id + '_T0'
ref.description = ref.description + \
    '. Added mutations accumulated before for the reference T0-WT.' 


## Write new sequence to file
with open(outfile, "w") as handle:
    SeqIO.write(ref, handle, "fasta")

print("File saved in: {}".format(outfile))


###########################################################################
##  TEST DUMMY EXAMPLES
###########################################################################

#### Sample string for working with indexes:

# asd = 'abcdefghijklmnopqrstuvwxyz'

# for i in range(0,len(asd)):
#     print('{} : {}'.format(i, asd[i]))

# ind=4       # INDEX = pos_str
# l=1         # LENGTH of the  pattern
# print(asd[ :ind] + 'E' + asd[ind + l: ])    #Slice before and after the required position

## With a longer pattern

# ind = 10
# l = 4
# asd[:ind] + 'KLMN' + asd[ind+l:]


#### Test sequence
# from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna
# my_dna = Seq("AGTACACTGGT", generic_dna)
