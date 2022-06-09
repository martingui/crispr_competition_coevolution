#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Alphabet import generic_dna 
from Bio.Blast.Applications import NcbiblastnCommandline
import pandas as pd
import sys
import os



# To launch this script in parallel:
# time find /home/enrique/work/Gandon/coevolution/bacter/Antoine/steps/blast_protospacers/samples/ -type d | parallel --bar --jobs 10 --max-proc 15 'python /home/enrique/work/Gandon/coevolution/bacter/Antoine/Script_pour_NGS_BIM/script_workflow/2_replace_protospacerv2_test_inter3.py {}'

######### FUNCTIONS

def do_blast(liste, blast_db):
    '''
    Gets a list of strings to be converted in sequences
    Writes all the sequences into a fasta file
    Makes blast outfmt6
    Loads results as pandas DF
    Gets only first match for every query
    Returns a list of subject names (sequence ID of the reference used)
    ''' 

    ### This drive was mounted to accelerate the writting and reading of fasta
    ### files and blast outputs. To make one use this commands in bash:
    # sudo mkdir -p /media/ramdisk
    # sudo mount -t tmpfs -o size=4096M tmpfs /media/ramdisk
    ### Else define another temporary directory
    
    temp_dir = '/media/ramdisk/' + dossier.replace('/home/enrique/work/Gandon/coevolution/bacter/Antoine/steps/blast_protospacers/samples/', '')

    fasta_file = temp_dir + '_'+ 'tmp.fasta'
    blast_out_file = temp_dir + '_'+ 'blast.out'
    
    # print(fasta_file, blast_out_file)

    create_fasta_sequences(liste, fasta_file)

    ## RUN BLAST
    blast_cline = NcbiblastnCommandline(query=fasta_file, db=blast_db, word_size=10, perc_identity=70, outfmt=6, out=blast_out_file)
    # blast_cline = NcbiblastnCommandline(query='tmp.fasta', db=blast_db, outfmt=6, num_alignments=1, out='blast.out')
    stdout, stderr = blast_cline()

    ## VARIABLE TO RETURN
    return_list = []

    ## IMPORT BLAST RESULTS
    ## IF RESULTS ARE EMPTY
    if os.stat(blast_out_file).st_size == 0 :
        ## Return -list with empty strings OR same list
        return liste
        # print("empty file")
        # pass
    else:
        ## GET BEST BLAST HIT (bbh) and return list of subject.id
        query_list, subject_list = bbh(blast_out_file)
        # return (query_list, subject_list)
    
    ## LOOK FOR ABSCENCE OF MATCHES IN LISTS WITH MULTIPLE ELEMENTS
        replacement_list = []
        for i in liste:
            if i in query_list:
                idx = query_list.index(i)
                subject_id_name = subject_list[idx].replace('_NC_007019.1_T0','')
                replacement_list.append(subject_id_name)
            else:
                replacement_list.append(i)
        # print(replacement_list)
        return replacement_list




def create_fasta_sequences(liste, fasta_out_file):
    '''
    Turns the sequences of one line into a fasta file
    '''
    ## CREATE FASTA FILES WHEREVER THIS RUNS
    records = []
    for i in range(0,len(liste)):
        # print(liste[0])
        record = SeqRecord(Seq(liste[i]), \
            description = 'seq_' + str(i), \
            id = liste[i] )
            # id = 'seq_' + str(i) )
        records.append(record)
    # print(record)

    with open(fasta_out_file, 'w') as handle:
        SeqIO.write(records, handle, 'fasta')
    

def bbh(blast_result):
    '''
    Get best blast hit
    '''

    blast_header = ['query_id', 'subject_id', 'identity', 'alignment_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bitScore']

    ## LOAD BLAST OUTPUT FILE
    with open(blast_result,'r') as handle:
        df = pd.read_csv(handle, sep="\t", names=blast_header)
        # print(df)

    ## TAKE ONLY THE FIRST LINE FOR EACH QUERY ID == TAKE ONLY BEST BLAST HIT
    df2 = df.groupby('query_id').first()
    # print(df2)

    ## SAVE TO FILE BBH
    # with open(inputFile+'.bbh', 'w') as outFile:
    #     df2.to_csv(outFile, sep='\t')

    return (list(df2.index), list(df2.subject_id))



#########

## PRINT WHEN IN SCRIPT




## GET ARGUMENTS AND PUT THEM IN VARIABLES

## DEACTIVATE BEFORE REAL RUN
# dossier = "/home/enrique/Documents/cefe/gandon/coevolution/bacter/Antoine/steps/antoine_scripts/B1T1"
# dossier = '/home/enrique/Documents/cefe/gandon/coevolution/bacter/Antoine/Script_pour_NGS_BIM/script_workflow'

## ACTIVATE BEFORE REAL RUN
# print('running...{0}  \n'.format(sys.argv[0]))
dossier = sys.argv[1]   # outdir
# dico = sys.argv[2]      # Fasta dictionary

## FILES

## eventually change paths
# os.chdir(path)
blast_db = "/home/enrique/work/Gandon/coevolution/bacter/Antoine/steps/blast_protospacers/blastdb/PAM_proto_blastdb"
input_to_replace = dossier + "/inter3"
# output_file = dossier + '/inter3_replaced'
output_file = dossier + "/inter4_replaced"
error_output = dossier + "/inter4_mismatch"

# print(dossier + '/'+ 'tmp.fasta')
# print(dossier + '/'+ 'blast.out')

# READ INPUT LINE BY LINE

print(dossier)

to_replace = []
with open(input_to_replace, 'r') as f:
    to_replace = f.read().splitlines()


## REMOVE SPACES ON THE LEFT. SPLIT INTO LIST
# for i in temp:
    # to_replace.append(i.lstrip().split(' '))


## PARSE THE LIST OF SEQUENCES

out_lines = []
for line in to_replace:
    # print(line)
    # counts = line[0]
    seq_list = line.replace('start-','').replace('-fin','')
    seq_list = seq_list.split('-')
    # print(seq_list)
    if seq_list[0] == 'fin':
        out_lines.append(line)
    else:
        # print(line)
        replacement_list = do_blast(seq_list, blast_db)
        # print(replacement_list)
        line = 'start-' + '-'.join(replacement_list) + '-fin'
        out_lines.append(line)

    # print(line)
    # print("#######################\n")

with open(output_file, 'w') as handle:
    for i in out_lines:
        handle.write(i + '\n')
