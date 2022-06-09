#! /home/enrique/envs/biopython/bin/python

## CONVERTS THE SEQUENCES FROM THE PROTOSPACERS TO FASTA FORMAT

### The file proto_sp_30pb_dico contains one sequence per line pf 30bp length  each
### It's their line number which is used as a reference for each

### IMPORT ALL THE LINES IN THE DICO FILE
handle = open('/mnt/alpha_raid/work/Gandon/phages/data/proto_sp_30pb_dico', 'r')
raw_sequences=handle.readlines()
handle.close()

### PARSE THE SEQUENCE'S LIST USING AN INDEX,
### USE THE INDEX OF EACH ELEMENT AS A NAME IN THE FASTA FORMAT
handle = open('/mnt/alpha_raid/work/Gandon/phages/data/proto_sp_30.fasta', 'w')
for i in range(1,len(raw_sequences) +1):
    handle.write('>proto-{}\n{}'.format(i, raw_sequences[i-1]))
handle.close()
