suppressPackageStartupMessages({
    library(data.table)
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

gene_seed_prob <- gene_seed_prob[gene_id %in% outliers$geneID]

fwrite(outliers, snakemake@output$outliers)
fwrite(gene_seed_prob, snakemake@output$gene_seed_prob_sub)
