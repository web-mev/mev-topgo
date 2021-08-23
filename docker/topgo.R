suppressMessages(suppressWarnings(library(org.Mm.eg.db)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
suppressMessages(suppressWarnings(library(topGO)))
suppressMessages(suppressWarnings(library(genefilter)))
suppressMessages(suppressWarnings(library(optparse)))

# args from command line:
args<-commandArgs(TRUE)

option_list <- list(
    make_option(
        c('-f','--input_file'),
        help='Path to the count matrix input.'
    ),
    make_option(
        c('-n', '--ontology'),
        help="The ontology class to analyze: MF, BP, or CC"
    ),
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
    ),
    make_option(
        c('-g', '--organism'),
        help='The organism database to use.'
    ),
    make_option(
        c('-s', '--identifier_type'),
        help='A string giving which tells us which identifier is used for the genes.'
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

# Check that the necessary inputs were provided:
if (is.null(opt$input_file)){
    message('Need to provide a count matrix with the -f/--input_file arg.')
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

if (! tolower(opt$identifier_type) %in% c("symbol", "ensembl", "entrez")) {
    message('For the gene identifier option, need to select one of: symbol, ensembl, or entrez.')
    quit(status=1)
}

# ingest the DGE output as a data frame
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
if(length(sigGenes) == 0){
    # Nothing to do since there were no significant genes at this threshold
    message('There were zero significant genes for the prescribed threshold of pvalue < ', opt$pvalue)

    # make dummy files and exit with zero since there isn't an error with the analysis
    # TODO: make files
    quit(status=0)
}


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

# an array of all gene IDs given by the original identifier
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
tryCatch({
        tgd <- new(
            "topGOdata",
            ontology = opt$ontology,
            allGenes = alg,
            nodeSize = opt$min_node_size,
            annot = annFUN.org,
            mapping = opt$organism,
            ID = opt$identifier_type
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
    }, error=function(x){
        message('Encountered an error when calculating GO enrichments. Often, this can be caused by specifying the incorrect gene identifier.')
        quit(status=1)
    }
)
# To the result table we'd like to add the actual 
# genes in each pathway. This way users can make
# genesets in the WebMeV UI.
mappings <- annFUN.org(opt$ontology, mapping = opt$organism, ID = opt$identifier_type)


# `mappings` is a list that looks like:
# $`GO:0000002`
#  [1] "ENSG00000151729" "ENSG00000025708" "ENSG00000068305" "ENSG00000115204"
#  [5] "ENSG00000198836" "ENSG00000196365" "ENSG00000117020" "ENSG00000275199"
#  [9] "ENSG00000114120" "ENSG00000140451" "ENSG00000171612" "ENSG00000125871"

# $`GO:0000003`
# [1] "ENSG00000147437" "ENSG00000125787" "ENSG00000189409" "ENSG00000183814"

# Turn that into a dataframe:
mapping_df = t(
    as.data.frame(
        lapply(mappings, function(x){
            return(paste(x, collapse=','))
        }),
        check.names=F
    )
)
colnames(mapping_df) <- c('genelist')

topgo.res = merge(topgo.res, mapping_df, by.x = 'GO.ID', by.y=0)


# Write the results to file
output_filename <- paste(
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
    row.names = F
)
