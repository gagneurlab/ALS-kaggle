suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})


outliers <- fread(snakemake@input$outliers)
gene_seed_prob <- fread(snakemake@input$gene_seed_prob)


outliers <- merge(
    outliers, 
    gene_seed_prob, 
    by.x = 'geneID',
    by.y = 'gene_id',
    all.x = TRUE
)

gene_seed_prob_sub <- gene_seed_prob[gene_id %in% outliers$geneID]

fwrite(outliers, snakemake@output$outliers)
fwrite(gene_seed_prob_sub, snakemake@output$gene_seed_prob_sub)

dt_gene_prob <- rbind(
    gene_seed_prob[, .(score, gene_set = 'genome')],
    gene_seed_prob_sub[, .(score, gene_set = 'expression outlier')]
)
dt_gene_prob[
    , gene_set := factor(
        gene_set,
        levels = c('genome', 'expression outlier'))
    ]

ggplot(dt_gene_prob, aes(-log10(score), fill = gene_set)) + 
    geom_histogram(position = 'identity', bins = 100, alpha = 0.8) +
    theme_classic() +
    theme(legend.position = 'bottom') +
    scale_fill_brewer(palette = 'Set1')

ggsave(snakemake@output$histogram)
