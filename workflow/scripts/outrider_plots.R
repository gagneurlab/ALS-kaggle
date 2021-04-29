suppressPackageStartupMessages({
    library(OUTRIDER)
    library(ggplot2)
})

ods <- readRDS(snakemake@input$ods)

png(snakemake@output$dim_search)
plotEncDimSearch(ods)
dev.off()

png(snakemake@output$power)
plotPowerAnalysis(ods)
dev.off()

png(snakemake@output$heatmap_raw)
plotCountCorHeatmap(ods, normalized = FALSE)
dev.off()

png(snakemake@output$heatmap_fit)
plotCountCorHeatmap(ods, normalized = TRUE)
dev.off()

png(snakemake@output$disp_est)
plotDispEsts(ods)
dev.off()

png(snakemake@output$per_sample)
plotAberrantPerSample(ods)
dev.off()

plotExpressedGenes(ods)
ggsave(snakemake@output$expr_genes)

res <- results(ods, all = TRUE)
ggplot(res, aes(pValue)) + geom_histogram(bins = 100)
ggsave(snakemake@output$pvalues)
