suppressMessages(suppressWarnings(library("org.Hs.eg.db", character.only=T, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library("org.Mm.eg.db", character.only=T, warn.conflicts = F, quietly = T)))

args<-commandArgs(TRUE)
organism <- args[1]
output_filename <- args[2]

if(organism == 'human'){
    db <- org.Hs.eg.db
    keys <- names(as.list(org.Hs.egALIAS2EG))
} else if(organism == 'mouse'){
    db <- org.Mm.eg.db
    keys <- names(as.list(org.Mm.egALIAS2EG))
} else {
    message('Unsupported organism choice.')
    quit(status=1)
}
# Create a mapping table. This will allow us to take the gene identifiers
# and map them to EntrezIDs for the fgsea process.
# ALIAS	ENSEMBL	SYMBOL	ENTREZID
# 1	A1B	ENSG00000121410	A1BG	1
# 2	A1B	ENSG00000172164	SNTB1	6641
# 3	ABG	ENSG00000121410	A1BG	1
# 4	GAB	ENSG00000121410	A1BG	1
# 5	HYST2477	ENSG00000121410	A1BG	1
gene_mappings <- select(
    db, 
    key=keys, 
    columns=c('ALIAS', 'ENSEMBL', 'SYMBOL', 'ENTREZID'), 
    keytype="ALIAS"
)
write.table(gene_mappings, output_filename, sep='\t', quote=T)
