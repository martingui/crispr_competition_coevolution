#!/bin/bash

#ligne a lancer dans terminal:
# 1- Commande; 2- dossier input; 3- dossier output; 4- fichier dico
# /home/manager/Documents/analyse_bioinfo/script_workflow/./run_all_files.sh       /home/manager/Documents/analyse_bioinfo/test_in/     /home/manager/Documents/analyse_bioinfo/test_out/    /home/manager/Documents/analyse_bioinfo/script_workflow/proto_sp_13pb_dico

 
#dossier contenant les données brutes. Doit finir par /
DIR=$1

#dossier de sortie. Doit finir par /
OUT=$2

#Dictionaire sous forme de liste (un par ligne) des spacers à remplacer
DICO=$3

pwd=$(dirname $0)

for i in $DIR*.fastq.gz
do 

    # echo $i;
    name=$(basename $i) ;
    echo 'Working on input file:' $name;
    shortname=${name%%_*};
    # echo $shortname;

    mkdir -p ${OUT}$shortname/;
    # echo ${OUT}$shortname/;


    # bash /home/manager/Documents/analyse_bioinfo/script_workflow/run_all_scripts.sh $i ${OUT}$shortname/ $DICO;
    bash $pwd/run_all_scripts.sh $i ${OUT}$shortname/ $DICO;

done


