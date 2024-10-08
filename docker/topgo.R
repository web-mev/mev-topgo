suppressMessages(suppressWarnings(library(org.Mm.eg.db)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
suppressMessages(suppressWarnings(library(topGO)))
suppressMessages(suppressWarnings(library(genefilter)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(jsonlite)))

# a function to handle writing of files. Used in a couple of places
write_results <- function(q, opt, working_dir) {
    # Write the results to file
    output_filename <- paste(
        "topGO",
        opt$ontology,
        "json",
        sep="."
    )
    output_filepath <- paste(working_dir, output_filename, sep='/')
    write(toJSON(q, auto_unbox=T), output_filepath)

    json_str = paste0(
        '{"go_results":"', output_filepath, '"}'
    )
    output_json <- paste(working_dir, 'outputs.json', sep='/')
    write(json_str, output_json)
} 

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
    message('Need to provide a feature table (e.g. a table of differentially expressed genes) with the -f/--input_file arg.')
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

if (! tolower(opt$identifier_type) %in% c("symbol", "ensembl", "refseq")) {
    message('For the gene identifier option, need to select one of: symbol, ensembl, or refseq.')
    quit(status=1)
}

# change the working directory to co-locate with the counts file:
working_dir <- dirname(opt$input_file)
setwd(working_dir)

# ingest the DGE output as a data frame
res <- read.table(
    file = opt$input_file,
    sep="\t",
    row.names = 1,
    header=T
)

# check that we have the proper columns:
if(!'padj' %in% names(res))
{
    message('The input file containing your significant genes (e.g. from a differential expression analysis) did not have a column named "padj", which we require for filtering.')
    quit(status=1)
}
if(!'overall_mean' %in% names(res))
{
    message('The input file containing your significant genes (e.g. from a differential expression analysis) did not have a column named "overall_mean", which we require. We use this column to create an appropriate "background" distribution of gene expressions for the topGO analysis.')
    quit(status=1)
}

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

    # make dummy data and exit with zero since there isn't an error with the analysis
    l=list(
        go_id='-',
        term='-',
        annotated=0,
        significant = 0,
        expected = 0,
        fisher_rank = 0,
        classic_pval = 0,
        elim_pval = 0,
        genelist=vector()
    )
    # needs to be a list in a list to export to the expected format.
    l = list(l)
    write_results(l, opt, working_dir)
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
        message('Encountered an error when calculating GO enrichments. Often, this can be caused by specifying the incorrect gene identifier. This can also be caused by choosing a threshold p-value that is too large and hence includes too many genes in the analysis.')
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


# `topgo.res` looks like:
# > head(topgo.res)
#        GO.ID                                 Term Annotated Significant
# 1 GO:0032501     multicellular organismal process      5909        3554
# 2 GO:0048731                   system development      3973        2440
# 3 GO:0099537             trans-synaptic signaling       615         458
# 4 GO:0007268       chemical synaptic transmission       608         453
# 5 GO:0098916 anterograde trans-synaptic signaling       608         453
# 6 GO:0099536                   synaptic signaling       633         468
#   Expected Rank in Fisher.classic Fisher.elim Fisher.classic
# 1  3219.62                      1     0.06047        < 1e-30
# 2  2164.76                      2     0.05637        1.0e-25
# 3   335.09                      3     0.28898        1.2e-25
# 4   331.28                      4     2.1e-05        1.8e-25
# 5   331.28                      5     1.00000        1.8e-25
# 6   344.90                      6     0.32988        4.9e-25

# Make the GO terms as the rownames
#rownames(topgo.res) <- topgo.res[,'GO.ID']

q = apply(topgo.res, 1, function(r){
    go_term = r[['GO.ID']]
    genelist = mappings[[go_term]]
    if(length(genelist) == 1){
        # add a NA (which gets translated to `null` by toJSON below)
        # so that the "unboxing" doesn't leave some genelists as
        # strings and others as lists of strings. By adding the null, GO
        # terms that correspond to single genes will look like:
        #{
        #     'GO:XXXXXX': {
        #         ...
        #         genelist: ['geneX', null]
        #     },
        #       ...
        # }
        genelist = c(genelist, NA)
    }
    # if the mappings list doesn't have any genes corresponding to a particular GO term,
    # then genelist is null
    if(is.null(genelist)){
        genelist = vector() # this way an empty array is encoded as "[]" in toJSON function
    }

    # Many of the terms are cast as strings-- change to what we expect
    annotated = as.integer(r[['Annotated']])
    significant = as.integer(r[['Significant']])
    fisher.rank = as.integer(r[['Rank in Fisher.classic']])
    expected = as.numeric(r[['Expected']])
    fisher.elim = as.numeric(r[['Fisher.elim']])
    # Note that p-values less than 1e-30 are given a p-value of "< 1e-30"
    # in the table. This causes parsing issues, so we cast to a hard zero here
    # The source code GenTable function above suggests we could pass a 'formatting function'
    # so that we can customize the printing, but this is a reasonable workaround. Directly
    # adding the `format.FUN` keyword argument caused errors
    if (trimws(r[['Fisher.classic']]) == '< 1e-30'){
        fisher.classic = 0
    } else {
        fisher.classic = as.numeric(r[['Fisher.classic']])
    }

    list(
        go_id=r[['GO.ID']],
        term=r[['Term']],
        annotated=annotated,
        significant = significant,
        expected = expected,
        fisher_rank = fisher.rank,
        classic_pval = fisher.classic,
        elim_pval = fisher.elim,
        genelist=genelist
    )
})

# this drops the names so that the json export does not have the GO terms as 'keys'.
# This makes it export as a list of objects
names(q) <- c()

write_results(q, opt, working_dir)
