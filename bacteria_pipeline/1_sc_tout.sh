#!/bin/bash

echo 'running...' $(basename $0) 

FILE=$1
outdir=$2

#suppression de la premiere ligne
zcat $FILE | sed 1d | 

#suppression de 3 ligne sur 4
sed -n '1~4p' > ${outdir}inter1_SeqTot
wc -l ${outdir}inter1_SeqTot > ${outdir}comptageinter1

#supprime ligne de plus de 590pb
awk 'length($0)<590' ${outdir}inter1_SeqTot > ${outdir}inter2_SeqSans590
wc -l ${outdir}inter2_SeqSans590 > ${outdir}comptageinter2

#Remplacement des repaet par un - en acceptant un mismatch
sed 's/[A-Z]TTGTACAGTTACTTAAATCTTGAGAGTACAAAAAC/-/g' ${outdir}inter2_SeqSans590 | 
sed 's/G[A-Z]TGTACAGTTACTTAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GT[A-Z]GTACAGTTACTTAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTT[A-Z]TACAGTTACTTAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTG[A-Z]ACAGTTACTTAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGT[A-Z]CAGTTACTTAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTA[A-Z]AGTTACTTAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTAC[A-Z]GTTACTTAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACA[A-Z]TTACTTAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAG[A-Z]TACTTAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGT[A-Z]ACTTAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTT[A-Z]CTTAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTA[A-Z]TTAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTAC[A-Z]TAAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACT[A-Z]AAATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTT[A-Z]AATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTA[A-Z]ATCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAA[A-Z]TCTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAA[A-Z]CTTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAAT[A-Z]TTGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATC[A-Z]TGAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCT[A-Z]GAGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTT[A-Z]AGAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTG[A-Z]GAGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTGA[A-Z]AGTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTGAG[A-Z]GTACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTGAGA[A-Z]TACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTGAGAG[A-Z]ACAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTGAGAGT[A-Z]CAAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTGAGAGTA[A-Z]AAAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTGAGAGTAC[A-Z]AAAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTGAGAGTACA[A-Z]AAAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTGAGAGTACAA[A-Z]AAC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTGAGAGTACAAA[A-Z]AC/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTGAGAGTACAAAA[A-Z]C/-/g' | 
sed 's/GTTGTACAGTTACTTAAATCTTGAGAGTACAAAAA[A-Z]/-/g' | 
#suppression des bases avant le premier -
sed 's/[^-]*\(-.*\)/\1/' | 
#Remplacement du premier spacers (commun a tous) par le mot start en acceptant un mismatch
sed 's/[A-Z]AATCTTGATTTGCTGTCAAACA/start/' |
sed 's/G[A-Z]ATCTTGATTTGCTGTCAAACA/start/' |
sed 's/GA[A-Z]TCTTGATTTGCTGTCAAACA/start/' |
sed 's/GAA[A-Z]CTTGATTTGCTGTCAAACA/start/' |
sed 's/GAAT[A-Z]TTGATTTGCTGTCAAACA/start/' |
sed 's/GAATC[A-Z]TGATTTGCTGTCAAACA/start/' |
sed 's/GAATCT[A-Z]GATTTGCTGTCAAACA/start/' |
sed 's/GAATCTT[A-Z]ATTTGCTGTCAAACA/start/' |
sed 's/GAATCTTG[A-Z]TTTGCTGTCAAACA/start/' |
sed 's/GAATCTTGA[A-Z]TTGCTGTCAAACA/start/' |
sed 's/GAATCTTGAT[A-Z]TGCTGTCAAACA/start/' |
sed 's/GAATCTTGATT[A-Z]GCTGTCAAACA/start/' |
sed 's/GAATCTTGATTT[A-Z]CTGTCAAACA/start/' |
sed 's/GAATCTTGATTTG[A-Z]TGTCAAACA/start/' |
sed 's/GAATCTTGATTTGC[A-Z]GTCAAACA/start/' |
sed 's/GAATCTTGATTTGCT[A-Z]TCAAACA/start/' |
sed 's/GAATCTTGATTTGCTG[A-Z]CAAACA/start/' |
sed 's/GAATCTTGATTTGCTGT[A-Z]AAACA/start/' |
sed 's/GAATCTTGATTTGCTGTC[A-Z]AACA/start/' |
sed 's/GAATCTTGATTTGCTGTCA[A-Z]ACA/start/' |
sed 's/GAATCTTGATTTGCTGTCAA[A-Z]CA/start/' |
sed 's/GAATCTTGATTTGCTGTCAAA[A-Z]A/start/' |
sed 's/GAATCTTGATTTGCTGTCAAAC[A-Z]/start/' |
#Remplacement de start et toutes les base le precedant par start
sed 's/.*start/start/' |
#Remplacement du debut du motif apres le dernier - par finseq
sed 's/CTCAAATGAA/finseq/' |
sed 's/[A-Z]TCAAATGAA/finseq/'|
sed 's/C[A-Z]CAAATGAA/finseq/'|
sed 's/CT[A-Z]AAATGAA/finseq/'|
sed 's/CTC[A-Z]AATGAA/finseq/'|
sed 's/CTCA[A-Z]ATGAA/finseq/'|
sed 's/CTCAA[A-Z]TGAA/finseq/'|
sed 's/CTCAAA[A-Z]GAA/finseq/'|
sed 's/CTCAAAT[A-Z]AA/finseq/'|
sed 's/CTCAAATG[A-Z]A/finseq/'|
sed 's/CTCAAATGA[A-Z]/finseq/'|
# Supression de seq et tout ce qui suit apres
sed 's/\(.*\)seq.*/\1/' |
#Trier en consevant tout ce qui commence par start
grep ^start |
#Trier en conservant tout ce qui se termine par fin
grep fin$ |
# Supprimer les lignes qui contiennent une suite de lettres spérieur a 32
sed -r '/[A-Z]{32,}/d' |
#supprimeles ligne qui contienne le mot finseq
sed -r '/finseq/d' |
# supprimer les lignes qui contiennent n'importe quelle lettre suivi de fin
sed -r '/[A-Z]fin/d' > ${outdir}inter3

wc -l ${outdir}inter3 > ${outdir}comptageinter3noreplaced

## Sorts alphabetically all the sequences in inter3
## counts how many times each is present
cat ${outdir}inter3 | sort | uniq -c > ${outdir}inter3sort
