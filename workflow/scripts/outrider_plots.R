suppressPackageStartupMessages({
    library(data.table)
    library(OUTRIDER)
    library(ggplot2)
})

anno_groups <- snakemake@params$anno_groups
ods <- readRDS(snakemake@input$ods)

png(snakemake@output$power)
plotPowerAnalysis(ods)
dev.off()

png(snakemake@output$per_sample)
plotAberrantPerSample(ods)
dev.off()

png(snakemake@output$qqplot)
plotQQ(ods, global = TRUE)
dev.off()

plotEncDimSearch(ods) + theme_bw()
ggsave(snakemake@output$dim_search)

plotExpressedGenes(ods)
ggsave(snakemake@output$expr_genes)

res <- results(ods, all = TRUE)
ggplot(res, aes(pValue)) + geom_histogram(bins = 100) + theme_classic()
ggsave(snakemake@output$pvalues)

cnts <- data.table(
    observed = as.vector(counts(ods)),
    expected = as.vector(normalizationFactors(ods)),
    aberrant = as.vector(aberrant(ods))
)
ggplot(cnts, aes(log10(expected + 1), log10(observed + 1), color = aberrant)) +
    geom_point() +
    scale_color_brewer(palette = 'Set1') +
    theme_classic() +
    theme(legend.position = 'bottom')
ggsave(snakemake@output$exp_obs_counts)

# cnts_long <- melt(
#     cnts,
#     measure.vars = c('expected', 'observed'),
#     variable.name = 'count_type',
#     value.name = 'counts'
# )
# 
# ggplot(cnts_long, aes(log10(counts + 1), after_stat(count), color = count_type)) +
#     geom_density(position = 'identity') +
#     scale_color_brewer(palette = 'Set1') +
#     theme_classic() +
#     theme(legend.position = 'bottom')
# 
