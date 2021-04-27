suppressPackageStartupMessages({
    library(OUTRIDER)
})

padj_cutoff <- snakemake@params$padj_cutoff

ods <- readRDS(snakemake@input$ods)
ods <- OUTRIDER(ods)
res <- results(ods, padjCutoff = padj_cutoff)

saveRDS(ods, snakemake@output$ods)
fwrite(res, snakemake@output$results)

