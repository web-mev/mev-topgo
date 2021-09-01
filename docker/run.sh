#!/bin/bash

DGE_FILE=$1
ORGANISM=$2
GENE_IDS=$3
ONTOLOGY=$4
PVALUE_THRESHOLD=$5
MIN_SIZE=$6
MAX_RESULT_SIZE=$7

if [ $ORGANISM = "Human" ]
then
    ORG_DB="org.Hs.eg.db"
elif [ $ORGANISM = "Mouse" ]
then
    ORG_DB="org.Mm.eg.db"
else
    echo "Not a valid organism choice." >&2
    exit 1;
fi

Rscript topgo.R \
    -f $DGE_FILE \
    -n $ONTOLOGY \
    -p $PVALUE_THRESHOLD \
    -g $ORG_DB \
    -s $GENE_IDS \
    -m $MIN_SIZE \
    -t $MAX_RESULT_SIZE