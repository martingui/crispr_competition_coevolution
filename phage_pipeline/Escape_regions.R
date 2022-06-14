
### This script makes the dataframe with start and end positions of escape region 
### An escape region is defined as a 37 nt long region, going from 1 nt before the PAM until the end of the protospacer

### Input:  PAM_protospacers.fasta
###           Also referred to as the PAM-protospacer dictionary

### Output: A dataframe saved as a csv, "PAM_positions.csv"
###           seq_name:     The name of PAM+protospacer as listed in the dictionary   
###           strand:       "-" or "+", referring to what strand the PAM is found on
###           sequence:     The 36 nt long PAM+protospacer sequence 
###           pos:          Position of the protospacer
###           escape_start: Start position of escape region 
###           escape_end:   End position of escape region


# Loading libraries
library("Biostrings")
library(tidyverse)
library(stringr)

## Loading the fasta file, aka the dictionary, and extracting features from it to be used in output dataframe
PAMS = readDNAStringSet("../data/PAM_protospacers.fasta")
header = names(PAMS)
seq_name = sub(" .*", "", header) 
strand <- str_match(header, "Strand=(.*?);") 
strand = strand[,2]
sequence = paste(PAMS)

## Making dataframe
df<-data.frame(seq_name, strand, sequence) 

## Extracting the genomic position of protospacer from seq_name and adding to "pos" column
pos <- df$seq_name %>% as.character() %>% str_match("PAM_(.*?)_NC") 
df$pos = pos[,2] %>% as.numeric()

## Finding start and end of escape region based on "pos" and sign of strand, i.e. "+" or "-".
## The region is 37 nt long 
df<- df %>% mutate(escape_start = ifelse(strand == "+", pos-32, pos-1))
df<- df %>% mutate(escape_end = ifelse(strand == "+", pos+4, pos+36))


## Save as csv file
write.csv(df, file = "../data/PAM_positions.csv", row.names = F)
