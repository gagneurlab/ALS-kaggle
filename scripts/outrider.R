library(OUTRIDER)

ctsTable <- read.table(snakemake@input[[1]], check.names=FALSE)

ods <- OutriderDataSet(countData=ctsTable)
ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
ods <- OUTRIDER(ods)
res <- results(ods)

fwrite(res, snakemake@output[[1]])

