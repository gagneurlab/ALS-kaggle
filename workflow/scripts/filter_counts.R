suppressPackageStartupMessages({
    library(ggplot2)
    library(OUTRIDER)
    library(GenomicFeatures)
})

fpkm_cutoff <- snakemake@params$fpkm_cutoff
#ensembl_release <- snakemake@params$ensembl_release

ctsTable <- read.table(
    snakemake@input$counts,
    header = TRUE,
    row.names = 1,
    sep = '\t'
)

# remove complete 0 rows
n_genes_unfiltered <- nrow(ctsTable)
ctsTable <- ctsTable[rowSums(ctsTable) > 0, ]
message(paste(
    n_genes_unfiltered - nrow(ctsTable),
    'removed that are not expressed in any sample')
)

# Create and filter Outrider dataset
ods <- OutriderDataSet(countData=ctsTable)
#txdb <- makeTxDbFromEnsembl(organism='Homo sapiens', release = ensembl_release)

ods <- computeGeneLength(ods, gtfFile = snakemake@input$gtf)
ods <- filterExpression(
    ods, 
    fpkmCutoff = fpkm_cutoff,
    gftFile = txdb,
    filterGenes = FALSE
)

png(snakemake@output$filtered_plot)
plotFPKM(ods)
dev.off()

saveRDS(ods, snakemake@output$ods)
