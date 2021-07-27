library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(topGO)
library(genefilter)
library(geneplotter)

# args from command line:
args<-commandArgs(TRUE)

option_list <- list(
    make_option(
        c('-f','--input_file'),
        help='Path to the count matrix input.'
    ),
    make_option(
        c('-o','--output_file_prefix'),
        help='The prefix for the output file'
    ),
    make_option(
        c('-n', '--ontology'),
        help="The ontology class to analyze: MF, BP, or CC"
    )
    make_option(
        c('-p','--pvalue'),
        default = 0.05,
        help='The p-value cutoff for significant genes'
    ),
    make_option(
        c('-m','--min_node_size'),
        default = 10,
        help='The minimum node size of a GO to consider for testing'
    ),
    make_option(
        c('-t', '--total_nodes'),
        default = 500,
        help='The total number of nodes to display in results'
    )
    make_option(
        c('-g', '--organism'),
        help='The organism database to use'
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

# Check that the necessary inputs were provided:
if (is.null(opt$input_file)){
    message('Need to provide a count matrix with the -f/--input_file arg.')
    quit(status=1)
}
if (is.null(opt$output_file_prefix)) {
    message('Need to provide the prefix for the output file with the -o arg.')
    quit(status=1)
}
if (is.null(opt$organism)) {
    message('Need to provide the organism with the -g/--organism arg.')
    quit(status=1)
}
if (is.null(opt$input_file)){
    message('Need to provide a ontology class with the -n/--ontology arg.')
    quit(status=1)
}

# Check if the ontology is a valid choice of BP, MF, or CC
if (! opt$ontology %in% c("MF", "BP", "CC")) {
    message('Need to provide a valid ontology of: MF, BP, or CC.')
    quit(status=1)
}


res <- read.table(
    file = opt$input_file,
    sep="\t",
    row.names = 1,
    header=T
)

# Get a list of the significant genes
sigGenes <- rownames(
    subset(
        res,
        padj < opt$pvalue
    )
)

# Get the mean values for all genes
overallBasemean <- as.matrix(res[, "overall_mean", drop = F])

# Get the indices of the significant genes in the mean matrix
sig_idx <- match(sigGenes, rownames(overallBasemean))

# Create the background list of genes matched by gene expression
backG <- c()
for (i in sig_idx){
    ind <- genefinder(overallBasemean, i, 10, method = "manhattan")[[1]]$indices
    backG <- c(backG, ind)
}
backG <- unique(backG)
# Convert list of genes to expressions
backG <- rownames(overallBasemean)[backG]

# Choice of ontologies
# This could be an input
onts = c( "MF", "BP", "CC" )

geneIDs = rownames(overallBasemean)
# Set universe for hypergeometric test as expression matched background
# and significant genes
inUniverse = geneIDs %in% c(sigGenes, backG)
# Input is the significant genes
inSelection = geneIDs %in% sigGenes
# Munging
alg <- factor(as.integer(inSelection[inUniverse]))
names(alg) <- geneIDs[inUniverse]

# prepare the data
tgd <- new(
    "topGOdata",
    ontology = opt$ontology,
    allGenes = alg,
    nodeSize = opt$min_node_size,
    annot = annFUN.org,
    mapping = opt$organism,
    ID = "symbol"
)
# Run tests
resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
resultTopGO.classic <- runTest(
    tgd, algorithm = "classic", statistic = "Fisher"
)
# append results
topgo.res <- GenTable(
    tgd,
    Fisher.elim = resultTopGO.elim,
    Fisher.classic = resultTopGO.classic,
    orderBy = "Fisher.classic",
    topNodes = opt$total_nodes
)

# Write the results to file
output_filename <- paste(
    opt$output_file_prefix,
    "topGO",
    opt$ontology,
    "tsv",
    sep="."
)
write.table(
    topgo.res,
    output_filename,
    sep="\t",
    quote=F,
    row.names = T
)
