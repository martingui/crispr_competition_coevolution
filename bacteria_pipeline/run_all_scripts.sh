#!/bin/bash

echo 'running... ' $(basename $0)


FILE=$1
outdir=$2
DICO=$3

pwd=$(dirname $0)

mkdir -p $outdir

## Line commented for update
bash $pwd/1_sc_tout.sh $FILE $outdir

python $pwd/2_replace_protospacerv2_inter3.py $outdir #$DICO

bash $pwd/2.5_sort_count_inter3.sh ${outdir}inter3_replaced $outdir





