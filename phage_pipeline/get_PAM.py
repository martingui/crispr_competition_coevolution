#!/usr/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna





#source file
# seq_test = "/home/manager/Documents/analyse_bioinfo/genome2972/seq_test"
# sf = "/home/manager/Documents/analyse_bioinfo/genome2972/Genome2972_1line.txt"
# fasta_file = 'data/refs/NC_007019.1_T0.fasta'
fasta_file='steps/new_ref_T0_cornellia_lotte/NC_007019.1_T0.fasta'
## This is the T0 file produced by Cornelia and Lotte which will be used from now on:
# fasta_file='people/lotte/coevolution/phages/data/refs/Sv2972'


#output destination file
# df = "/home/manager/Documents/analyse_bioinfo/genome2972/output.txt"
df = 'steps/PAM_detection/PAM_protospacers.fasta'


## MOTIF 
pat = "AGAA"
pat2 = "GGAA"

## Motif as sequence object
pat_seq = Seq(pat, generic_dna)
pat2_seq = Seq(pat2, generic_dna)
## Protospacer length is 30 + 2 nuc betweent the PAM and the protospacer
prot_len = 32

## FASTA HEADER INFORMATION
## Fasta ID will be: PAM_CoordRef_genome
## CoordRef is the index of the left most letter of the motif
id_pre = 'PAM'
id_post = fasta_file.replace('steps/new_ref_T0_cornellia_lotte/','').replace('.fasta', '')

## The description for the fasta record:
## Motif=xxxx;Strand=+;Protospacer_length=N
desc_motif = 'Motif='
desc_strand = 'Strand='
desc_proto_len = 'Protospacer_lenth=' + str(prot_len)


## Import sequence
# with open(sf, 'r') as f:
#     texte = f.readlines()[0]
#     texte = texte.upper()

ref = SeqIO.read(fasta_file, "fasta")


## SLIDING WINDOW
indices = []
records = []


count_pat_for=0     # patter1 forward
count_pat_crv=0     # pattern1 reverse complementary
count_pat2_for=0    # pattern2 forward
count_pat2_crv=0    # pattern2 reverse complementary

for i in range(0,len(ref.seq)):
    a = ref.seq[i : i + len(pat)]

    if pat_seq in a:
        count_pat_for +=1
        indices.append(i)
        prot_PAM_seq = ref.seq[i - prot_len : i + len(pat_seq)]
        record = SeqRecord( seq=prot_PAM_seq, \
            id = '_'.join([id_pre, str(i+1), id_post]), \
            description =  ';'.join([desc_motif+pat, desc_strand+'+', desc_proto_len]) )
        records.append(record)

    elif pat_seq.reverse_complement() in a:
        count_pat_crv +=1
        indices.append(i)
        prot_PAM_seq = ref.seq[i : i + len(pat_seq) + prot_len]
        record = SeqRecord(seq=prot_PAM_seq, \
            id = '_'.join([id_pre, str(i+1), id_post]), \
            description =  ';'.join([desc_motif+pat, desc_strand+'-', desc_proto_len]))
        records.append(record)

    elif pat2_seq in a:
        count_pat2_for +=1
        indices.append(i)
        prot_PAM_seq = ref.seq[i - prot_len : i + len(pat2_seq)]
        record = SeqRecord( seq=prot_PAM_seq, \
            id = '_'.join([id_pre, str(i+1), id_post]), \
            description =  ';'.join([desc_motif+pat2, desc_strand+'+', desc_proto_len]) )
        records.append(record)

    elif pat2_seq.reverse_complement() in a :
        count_pat2_crv +=1
        indices.append(i)
        prot_PAM_seq = ref.seq[i : i + len(pat2_seq) + prot_len]
        record = SeqRecord(seq=prot_PAM_seq, \
            id = '_'.join([id_pre, str(i+1), id_post]), \
            description =  ';'.join([desc_motif+pat2, desc_strand+'-', desc_proto_len]))
        records.append(record)


with open(df, 'w') as handle:
    SeqIO.write(records, handle, 'fasta')

print("The Motif {0} was found {1} times in the forward strand \
and {2} times in the reverse strand".format(pat, \
count_pat_for, \
count_pat_crv))


print("The Motif {0} was found {1} times in the forward strand \
and {2} times in the reverse strand\n".format(pat2, \
count_pat2_for, \
count_pat2_crv))

print("There are in total {} fasta sequences".format(len(records)))
print("The output file is: {}".format(df))

# with open(df, 'w') as dest:
#     for i in indices:
#         if i < 32: 
#             pass
#         else:
#             output = texte[i - 32 : i + len(pat)]
#             print(output)
#             dest.write(output + "\n")
