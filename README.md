# mev-topgo

This repository contains a WebMeV-compatible Docker image and scripts for executing a gene ontology analysis as performed by Bioconductor's topGO package (https://bioconductor.org/packages/release/bioc/html/topGO.html)

This implementation assumes you have previously run a differential gene expression analysis and exported a file that indicates genes which are differentially expressed for a provided contrast (e.g. wild-type versus knockout)

Inputs to the script include:
- A table (tab-delimited format) providing information about differential expression. To enable seamless integration of different WebMeV tools, we impose some restrictions on this file which are consistent with the format of the outputs from our differential expression testing tools:
    - Must include a column named "overall_mean" which can simply be the mean across all samples, regardless of their experimental group. The idea is that we use this column to create a "background" distribution of similarly expressed genes. That background/null distribution will then be used to test for enrichment of genes related to a particular GO-term.
    - Must include a "padj" column, which can either be the adjusted p-value or just an unadjusted p-value column.

Other, more basic options are provided below.

---

### To run external of WebMeV

Either:
- build the Docker image using the contents of the `docker/` folder (e.g. `docker build -t myuser/topgo:v1 .`) 
- pull the docker image from the GitHub container repository (see https://github.com/web-mev/mev-topgo/pkgs/container/mev-topgo)

To run, change to the directory where your differential expression results are located. Then:
```
docker run -it -v $PWD:/work Rscript /usr/local/bin/topgo.R \
    -f /work/<differential expression results> \
    -n <The GO terms to test. One of MF,CC,BP> \
    -p <threshold for removing genes with values of padj greater than this value> \
    -g <organism. Either "org.Hs.eg.db" or "org.Mm.eg.db"> \
    -s <gene identifier type. One of symbol, ensembl, or refseq> \
    -m <minimum size for GO term>\
    -t <Max number of results to output>
```
Some notes:
- MF, CC, BP control which ontology we look at and test. 
    - MF: molecular function
    - CC: cellular component
    - BP: biological process
- We only allow GO enrichment testing for human and mouse, which have the most well-annotated databases. The "org.Hs.eg.db" refers to human, and "org.Mm.eg.db" is for mouse.
- To look up the gene symbol and potentially map to the database identifiers, we need to know which annotation system you are using:
    - "symbol" refers to the "common" gene name (e.g. TP53, KRAS, etc.)
    - "ensmebl" refers to Ensembl-based identifiers, which start with ENSG for human and ENSMUSG for mouse
    - "RefSeq" are identifiers that start with "NM", "NP", "XM", etc.