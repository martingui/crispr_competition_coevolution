## Description of the workflow.

### Run the workflow

To run this workflow on multiple fastq.gz files run the following command:
Replace the variables according to your system.

```bash
inputdir=/home/enrique/work/Gandon/coevolution/bacter/Antoine/Script_pour_NGS_BIM/test_in/
outputdir=/home/enrique/work/Gandon/coevolution/bacter/Antoine/steps/antoine_scripts
protospacer_file=/home/enrique/work/Gandon/coevolution/bacter/Antoine/steps/PAM_detection/PAM_protospacers.fasta

bash run_all_files.sh $inputdir  $outputdir  $protospacer_file
```

----


### Description

What does each script ? -- The aim of this workflow was to identify the spacers of the bacteria

The data is long single end reads which cover the whole studied region.
There are a lot of text replacements realized by the script `1_sc_tout.sh`,
this allows a manual control of the results and makes it easier to have a grip on what we want to count.


### Run in all input files

The script `run_all_files.sh` makes a loop on all the input **fastq.gz** and launches `run_all_scripts.sh`. 
It also creates an output directory per sample, rename your input files to reflect biological siginificant names.

The script `run_all_scripts.sh` launches each of the other scripts by passing the output directory.

### `1_sc_tout.sh` 

selects only the sequence of the fastq.gz file and does a series of string replacements:

1. Replaces the *repeat* sequence with a *-*. Pattern: **GTTGTACAGTTACTTAAATCTTGAGAGTACAAAAAC** max 1 mismatch.
2. The bases before the *repeat* are removed.
3. The first *spacer* is replaced by *start*. Pattern: **GAATCTTGATTTGCTGTCAAACA** max 1 mismatch.
4. Replaces all the *start* and all the preceding bases by *start*, allows to visually align the sequences
5. Replaces the bases after the last repeat (*-*) by *finseq*. Pattern **CTCAAATGAA** max 1 mismatch.
6. Removes *seq* until the end of the string.
7. Selects only the strings starting with *start* and ending with *fin*
8. Removes sequences longer than 32
9. Removes the sequences containing the pattern *finseq* -- There should be none.
10. Removes sequences any sequences containing any letter followed by *fin*.
11. Creates an intermediary output file: **inter3**


### `2_replace_protospacerv2_inter3.py`

This scripts blasts each spacer (query) and replaces it with the best blast hit if there is a match. Runs slow as each sequence is blasted, a better version has been developped for later uses.


### `2.5_sort_count_inter3.sh`

This script sort and counts the lines of inter3 to regroup the same spacers and get the count of their occurrences 



### `Bact_genos.py`
Reads in bacterial genotypes and uses their Id's to locate the middle position of the escape region each spacer targets. 
Makes three outputs: 
  + data/bateria/Bacteria_genos.csv
  + data/bacteria/no_match_R.txt
  + data/bacteria/no_match_W.txt



