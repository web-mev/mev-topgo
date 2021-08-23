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
        c('-o','--output_file_prefix'),
        help='The prefix for the output file'
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
        c('-g', '--gene_mapping'),
        help='A file giving the various gene symbols and mapping between them'
    ),
    make_option(
        c('-j', '--go_mapping'),
        help='A file giving the GO terms mapped to the associated genes.'
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

if (is.null(opt$output_file_prefix)) {
    message('Need to provide the prefix for the output file with the -o arg.')
    quit(status=1)
}

if (is.null(opt$input_file)){
    message('Need to provide a ontology class with the -n/--ontology arg.')
    quit(status=1)
}

if (is.null(opt$gene_mapping)){
    message('Need to provide a file with gene symbol mapping via the -g/--gene_mapping arg.')
    quit(status=1)
}

if (is.null(opt$go_mapping)){
    message('Need to provide a file mapping GO terms to genes via the -j/--go_mapping arg.')
    quit(status=1)
}

# Check if the ontology is a valid choice of BP, MF, or CC
if (! opt$ontology %in% c("MF", "BP", "CC")) {
    message('Need to provide a valid ontology of: MF, BP, or CC.')
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

print(head(sigGenes))

# Get the mean values for all genes
overallBasemean <- as.matrix(res[, "overall_mean", drop = F])
print(head(overallBasemean))
print(dim(overallBasemean))

# Get the indices of the significant genes in the mean matrix
sig_idx <- match(sigGenes, rownames(overallBasemean))

# Create the background list of genes matched by gene expression
backG <- c()
for (i in sig_idx){
    ind <- genefinder(overallBasemean, i, 10, method = "manhattan")[[1]]$indices
    backG <- c(backG, ind)
}
backG <- unique(backG)
print(head(backG))
# Convert list of genes to expressions
backG <- rownames(overallBasemean)[backG]
print(head(backG))
print('????')

# an array of all gene IDs given by the original identifier
geneIDs = rownames(overallBasemean)

# Now we map to Ensembl IDs
GENE_ID_TYPE <- tolower(opt$identifier_type)
if (GENE_ID_TYPE != 'ensembl') {
    if(GENE_ID_TYPE == 'symbol'){
        # note that when we map back at the end, we use the symbol column
        # since the many:1 mapping of alias to entrez causes problems.    
        chosen_col = 'ALIAS'
        remap_col = 'SYMBOL'
    } else if(GENE_ID_TYPE == 'entrez'){
        chosen_col = 'ENTREZID'
        remap_col = 'ENTREZID'
    } else {
        message('Could not understand which gene identifier to map from.')
        quit(status=1)
    }

    # merge to keep only those where we have a mapping:
    gene_info_df = read.table(opt$gene_mapping)
    gene_info_df = unique(gene_info_df[gene_info_df[,chosen_col] %in% geneIDs, c(chosen_col,'ENSEMBL', remap_col)])

    # If no remaining rows, error out
    if(dim(gene_info_df)[1] == 0){
        message('After mapping the gene identifiers, there were no remaining rows. Was the choice of gene identifier correct?')
        quit(status=1)
    }
    # Note that a single gene symbol can map to multiple Ensembl IDs. Such as this:
    #         ALIAS         ENSEMBL
    # 1552    APLP2 ENSG00000084234
    # 1950    ASAH1 ENSG00000104763
    # 4357     CD28 ENSG00000178562
    # 7093     CTSS ENSG00000163131
    # 8853     TYMP ENSG00000025708
    # 8854     TYMP ENSG00000284194
    # 12789     GRN ENSG00000030582
    # 24615     CFP ENSG00000126759
    # 26757    PSAP ENSG00000197746
    # 26758    PSAP ENSG00000137409
    # 26759    PSAP ENSG00000122852
    # 26760    PSAP ENSG00000185303
    # 33471   TGFBI ENSG00000120708
    # 49040   RCAN3 ENSG00000117602
    # 53362 ZNF385A ENSG00000161642
    # 69446  SCPEP1 ENSG00000121064
    # 74047   TTYH3 ENSG00000136295
    # 97213 FAM102A ENSG00000167106
    backG <- gene_info_df[gene_info_df[,chosen_col] %in% backG, 'ENSEMBL']
    sigGenes <- gene_info_df[gene_info_df[,chosen_col] %in% sigGenes, 'ENSEMBL']
    geneIDs <- gene_info_df[,'ENSEMBL']
} 
print(head(backG))
print('------------------------')
print(head(sigGenes))

# Set universe for hypergeometric test as expression matched background
# and significant genes
inUniverse = geneIDs %in% c(sigGenes, backG)
# Input is the significant genes
inSelection = geneIDs %in% sigGenes
# Munging
alg <- factor(as.integer(inSelection[inUniverse]))
names(alg) <- geneIDs[inUniverse]
print(head(alg))

# read the proper mapping of GO terms to Ensembl genes:
go2genes <- readMappings(opt$go_mapping)

print(head(go2genes))
# prepare the data
tgd <- new(
    "topGOdata",
    ontology = opt$ontology,
    allGenes = alg,
    nodeSize = opt$min_node_size,
    annot = annFUN.GO2genes,
    GO2genes=go2genes
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
print(head(topgo.res))

# To the result table we'd like to add the actual 
# genes in each pathway. This way users can make
# genesets in the WebMeV UI.
go_terms = topgo.res[,'GO.ID']
go2genes_subset = go2genes[go_terms] # gg

# if the original table did not have Ensembl IDs,
# map to that original gene identifier system
if (GENE_ID_TYPE != 'ensembl') {

}
go_genes_df = t(
    as.data.frame(
        lapply(go2genes[go_terms], function(x){
            return(paste(x, collapse=','))
        }),
        check.names=F
    )
)
colnames(go_genes_df) = c('genelist')

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
