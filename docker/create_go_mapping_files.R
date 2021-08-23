suppressMessages(suppressWarnings(library("org.Hs.eg.db", character.only=T, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library("org.Mm.eg.db", character.only=T, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(topGO)))

args<-commandArgs(TRUE)
organism <- args[1]
ontology <- args[2]
output_filename <- args[3]

if(organism == 'human'){
    db <- 'org.Hs.eg.db'
} else if(organism == 'mouse'){
    db <- 'org.Mm.eg.db'
} else {
    message('Unsupported organism choice.')
    quit(status=1)
}

# Use the Ensembl annotations
mappings <- annFUN.org(ontology, mapping = db, ID = "Ensembl")


# `mappings` is a list that looks like:
# $`GO:0000002`
#  [1] "ENSG00000151729" "ENSG00000025708" "ENSG00000068305" "ENSG00000115204"
#  [5] "ENSG00000198836" "ENSG00000196365" "ENSG00000117020" "ENSG00000275199"
#  [9] "ENSG00000114120" "ENSG00000140451" "ENSG00000171612" "ENSG00000125871"

# $`GO:0000003`
# [1] "ENSG00000147437" "ENSG00000125787" "ENSG00000189409" "ENSG00000183814"

# Turn that into a dataframe:
df = t(
    as.data.frame(
        lapply(mappings, function(x){
            return(paste(x, collapse=','))
        }),
        check.names=F
    )
)
write.table(df, output_filename, sep='\t', quote=F, col.names=F)