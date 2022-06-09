#! /bin/bash




blastn -task blastn-short \
	-db /mnt/alpha_raid/work/Gandon/phages/refs/Sv2972/NC_007019.1.fasta \
	-query /mnt/alpha_raid/work/Gandon/phages/data/proto_sp_13.fasta \
	-out /mnt/alpha_raid/work/Gandon/phages/results/blast_protospacers/blast_proto13_Sv.tab \
	-outfmt 7 -num_threads 30 \
	-word_size 13

## NOTES
# There are 694 protospaces of 13pb each
# On the first try I didn't used the -word_size parameter (default7)
# 	There were too many results, with mostly only one with a 13pb matchs and 100%
# To diminish the results I added the parameter -word_size 13.
# This outputs 689 protospacers. 5 are not matching.
# 	-- > I'll try to get them on another blast, 

# 	FOR NOW i'll try to do the freebayes part to put all the clones and T0 in only one VCF with Freebayes


blastn -task blastn-short \
	-db /mnt/alpha_raid/work/Gandon/phages/refs/Sv2972/NC_007019.1.fasta \
	-query /mnt/alpha_raid/work/Gandon/phages/data/mix_BIM_17locus.fasta \
	-out /mnt/alpha_raid/work/Gandon/phages/results/blast_protospacers/blast_mix_BIM_17_loci.tab \
	-outfmt 7 -num_threads 30 \
	-word_size 13


blastn -task blastn-short \
	-db /mnt/alpha_raid/work/Gandon/phages/data/mix_BIM_17locus.fasta \
	-query /mnt/alpha_raid/work/Gandon/phages/data/proto_sp_13.fasta \
	-out /mnt/alpha_raid/work/Gandon/phages/results/blast_protospacers/blast_protospacer-mix_BIM_17_loci.tab \
	-outfmt 7 -num_threads 30 \
	-word_size 12



blastn -task blastn-short \
	-db /mnt/alpha_raid/work/Gandon/phages/refs/Sv2972/NC_007019.1.fasta \
	-query /mnt/alpha_raid/work/Gandon/phages/data/mix_BIM_17locus.fasta \
	-out /mnt/alpha_raid/work/Gandon/phages/results/blast_protospacers/blast_mix_BIM_17_loci.tab \
	-outfmt 7 -num_threads 30 \
	-word_size 12


blastn -task blastn-short \
	-db /mnt/alpha_raid/work/Gandon/phages/refs/Sv2972/NC_007019.1.fasta \
	-query /mnt/alpha_raid/work/Gandon/phages/data/proto_sp_13.fasta \
	-out /mnt/alpha_raid/work/Gandon/phages/results/blast_protospacers/test.tab \
	-outfmt 7 -num_threads 30 \
	-word_size 12


blastn -task blastn-short \
	-db /mnt/alpha_raid/work/Gandon/phages/refs/Sv2972/NC_007019.1.fasta \
	-query /mnt/alpha_raid/work/Gandon/phages/data/proto_sp_30.fasta \
	-out /mnt/alpha_raid/work/Gandon/phages/results/blast_protospacers_30/blast_proto30-genome.tab \
	-outfmt 6 -num_threads 30 \
	-word_size 20


blastn -task blastn-short \
	-db /mnt/alpha_raid/work/Gandon/phages/data/mix_BIM_17locus.fasta \
	-query /mnt/alpha_raid/work/Gandon/phages/data/proto_sp_30.fasta \
	-out /mnt/alpha_raid/work/Gandon/phages/results/blast_protospacers_30/blast_proto30-mix_BIM.tab \
	-outfmt 6 -num_threads 30 \
	-word_size 28

