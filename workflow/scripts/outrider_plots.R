suppressPackageStartupMessages({
    library(data.table)
    library(OUTRIDER)
    library(ggplot2)
    library(viridis)
    library(patchwork)
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

obs_cnt_mtx <- counts(ods, normalized = FALSE)
exp_cnt_mtx <- normalizationFactors(ods)

obs_cnt_mtx_log <- log10(obs_cnt_mtx + 1)
exp_cnt_mtx_log <- log10(exp_cnt_mtx + 1)

cnts <- data.table(
    # geneID=rep(rownames(ods), ncol(ods)), 
    # sampleID=rep(colnames(ods), each=nrow(ods)),
    observed = as.vector(obs_cnt_mtx_log) - rowMeans(obs_cnt_mtx_log),
    expected = as.vector(exp_cnt_mtx_log) - rowMeans(exp_cnt_mtx_log),
    aberrant = as.vector(aberrant(ods))
)

# r_squared <- round(summary(lm(observed ~ expected, data = cnts))$r.squared, 2)

p1 <- ggplot(cnts, aes(x=expected, y=observed)) +
    geom_hex(
        aes(fill=stat(log10(count))),
        binwidth = 0.02,
        show.legend = FALSE
    ) +
    scale_fill_gradientn(colours = viridis(30)) +
    geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") +
    # annotate(
    #     'text',
    #     x = -1,
    #     y = 2,
    #     label = paste0("r^2 == ", r_squared),
    #     parse = TRUE,
    #     size = 4
    # ) +
    labs(x="expected log10 counts", y="observed log10 counts") +
    theme_classic()


obs_cnt_mean_mtx <- log10(rowMeans(obs_cnt_mtx) + 1)
exp_cnt_mean_mtx <- log10(rowMeans(exp_cnt_mtx) + 1)
cnts_mean <- data.table(
    observed = obs_cnt_mean_mtx - mean(obs_cnt_mean_mtx),
    expected = exp_cnt_mean_mtx - mean(exp_cnt_mean_mtx)
)

# p2 <- ggplot(cnts_mean, aes(x=expected, y=observed)) +
#     geom_hex(
#         aes(fill=stat(log10(count))),
#         binwidth = 0.1,
#         show.legend = FALSE
#     ) +
#     scale_fill_gradientn(colours = viridis(30)) +
#     geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") +
#     labs(x="expected log10 counts", y="observed log10 counts") +
#     theme_classic()
# 
# p1 + p2
p1
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
