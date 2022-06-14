## Contents of scripts folder


Files:

* 00_create_py_env.sh
* 01_quality_check.sh*
* 02_trimm_and_clean.sh*
* 03_mapping.sh*
* 04_snpcalling.sh*
* 05b_convert_protospacer_dico2fasta.py*
* 06b_blast_protospaces.sh*

* procedure.sh
* README.md
* requirements_py-env.txt

* vcf_parser3.py
* vcf_to_csv.py
* ReadData_Functions.R
* Make_VCM_data.R


----


## File descriptions


### 00_create_py_env.sh

Creates a python virtual environment using `virtualenv`, the default python3 version of the system and will storte the environment in `~/envs/coev`. The installation of packages is done through pip.


### 01_quality_check.sh*

It will use FastQC to create quality control reports and then use multiqc to assemble the reports in only file. To make things easier, the input files are separated in 3 groups R, W and Other. These groups come from different treatments.

This script takes one argument: The path to the working directory, which is the project directory: `/home/user/work/coevolution/phages/`. **Don't forget that final stroke** 


### 02_trimm_and_clean.sh*

Launches Trimmomatic to clean data.
The parameters are embeded in the code

This script takes one argument: The path to the working directory, which is the project directory: `/home/user/work/coevolution/phages/`. **Don't forget that final stroke** 


### 03_mapping.sh*

The index creation is commented in the top of the script. It's only required once.
The path to the input files is a full path using a variable.
The mapper is bowtie2, after mapping the sam is sorted and converted to a bam and indexed so it's ready for the next stage.

This script takes one argument: The path to the working directory, which is the project directory: `/home/user/work/coevolution/phages/`. **Don't forget that final stroke** 


### 04_snpcalling.sh*

It uses freebayes to make the snp calling

This script takes three arguments: 
1 The path to the working directory, which is the project directory: `/home/user/work/coevolution/phages`. **Don't forget that final stroke** 
2 Path to the reference
3 Path to the output directory


### 05b_convert_protospacer_dico2fasta.py*

The *protospacer_dico* is a raw text file with only one sequence per line.
Each line is a protospacer sequence that has been manually selected and curated.
This script converts these sequences to fasta format.
The name given to each sequence is the line number.

No arguments. the path to the file is hard written


### 06b_blast_protospaces.sh*

Protospacers come from bacteria, we need to know where do they map on the virus' genome.

We tried out different parameters, this script is more a reminder of the commands and parameters tested.



----



## The other scripts

### procedure.sh

The commands used to launch the scripts up here as well as the supplementary commands to separate the different data, extract and all other action is written here.

It contains a pre-treatment of data to create sub-groups using symbolic links.



### update_ref_genome.py

Allows to create a consensus fasta file from the reference fasta
and a VCF file from which the variants will be integrated.
It filters the variants with a minimum frequency of 0.45.

It is equivalent to bcftools consensus:
`bcftools consensus --sample unknown -f NC_007019.1.fasta TO-WT_S83.vcf.gz -o test.fasta`


### get_PAM.py -- Getting new PAMs

From the genome of Streptococcus virus 2972 we'll be looking for new PAMs (Protospacer Associated Motif).

The researched patterns are: **AGAA** and **GGAA**; as well as the reverse complementary sequences **TTCT** and **TTCC**.

The PAMS are little sequences of 4 nucleotides located on the 3' side of a protospacer sequence. The length of a protospacer taken into account is 32 + 4 (protospacer + PAM).


#### Example of sequences & Nomenclature

Here are the templates of sequences.

The PAM sequences are found in 3 prime of the protospacers.
Each needs a unique name of ID for the fasta header.

```
>PAM_Coord_ref_genome Motif=AGAA;Strand=+;Protospacer_length=N
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxAGAA
>PAM_Coord_ref_genome Motif=GGAA;Strand=+;Protospacer_length=N
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxGGAA
>PAM_Coord_ref_genome Motif=AGAA;Strand=-;Protospacer_length=N
TTCTxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
>PAM_Coord_ref_genome Motif=GGAA;Strand=-;Protospacer_length=N
TTCCxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

The nomenclature of the fasta header goes as follows:

```
>PAM_152_Sv2972 M=AGAA;Strand=+;Protospacer_length=32
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxAGAA
>PAM_10863_Sv2972 M=TTCT;Strand=-;Protospacer_length=32
TTCTxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

Mandatory fields for ID field in header FASTA. No spaces

* **PAM** = Type of sequence **P**rotospacer **A**ssociated **M**otif
* **152** = The coordinate of the first nucleotide of the motif.
* **Sv2972** = Reference genome used to find these sequences

Comment field for FASTA format. No strict rules. I used Definition=Value without spaces. Each separated by a semi-colon

* **Motif=AGAA** = Motif The motif used to find that sequence
* **Strand=+** = The strand where the sequence is. If the motif is on the strand -, the sequence present will be the reverse complementary.
* **Protospacer_length** = The lenght of sequence concidered to be a protospacer without the lenght of the PAM.


### Escape_regions.R

Uses the output of get_PAM.py to produce a file with the positions of all PAMS: PAM_positions.csv


### requirements_py-env.txt

File used by the script `00_create_py_env.sh`
It contains the list of libraries and versions used in python for this analysis.


### vcf_parser3.py

It contains allows to import data from vcf files and put it in a pandas dataframe format to ease the analyses.

It works with objects: The _attributes_ contain information and the _methods_ are internal functions for that object. For example "remove control SNPs"

The last class allows to create a heatmap with matplotlib.

### vcf_to_csv.py

Function that takes the .vcf files from FreeBayes and outputs a .csv
  
### ReadData_Functions.R**

Functions that read and tidy up the csv files produced by "vcf_to_csv.py". 

### Make_VCM_data.R

Function that makes .tsv data for the analysis from the .csv created in 'vcf_to_csv.py'
Also executes the function and makes the data.
  




