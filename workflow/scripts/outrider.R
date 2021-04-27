suppressPackageStartupMessages({
    library(OUTRIDER)
})

padj_cutoff <- snakemake@params$padj_cutoff

ods <- readRDS(snakemake@input$ods)

# apply filter
ods <- ods[mcols(ods)$passedFilter,] 

# run end-to-end pipeline
ods <- OUTRIDER(ods)
saveRDS(ods, snakemake@output$ods)

# extract and save results
res <- results(ods, padjCutoff = padj_cutoff)
fwrite(res, snakemake@output$results)

