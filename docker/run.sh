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

if [ "$ONTOLOGY" = "MF (Molecular Function)" ]
then
    ONT_SHORT="MF"
elif [ "$ONTOLOGY" = "CC (Cellular Component)" ]
then
    ONT_SHORT="CC"
elif [ "$ONTOLOGY" = "BP (Biological Process)" ]
then
    ONT_SHORT="BP"
else
    echo "Not a valid ontology choice." >&2
    exit 1;
fi

Rscript /usr/local/bin/topgo.R \
    -f $DGE_FILE \
    -n $ONT_SHORT \
    -p $PVALUE_THRESHOLD \
    -g $ORG_DB \
    -s $GENE_IDS \
    -m $MIN_SIZE \
    -t $MAX_RESULT_SIZE
