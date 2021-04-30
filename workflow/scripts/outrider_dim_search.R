suppressPackageStartupMessages({
    library(data.table)
    library(OUTRIDER)
})

total_threads <- snakemake@threads

ods <- readRDS(snakemake@input$ods)

# apply filter
ods <- ods[mcols(ods)$passedFilter,]

# encoding dimension search
findEncodingDim_opt <- function (
    ods,
    params,
    freq = 0.01,
    zScore = 3,
    sdlog = log(1.6),
    lnorm = TRUE,
    inj = "both",
    ...,
    BPPARAM = bpparam()
) {
    ods <- estimateSizeFactors(ods)
    ods <- OUTRIDER:::injectOutliers(
        ods,
        freq = freq,
        zScore = zScore,
        inj = inj,
        lnorm = lnorm,
        sdlog = sdlog
    )
    eval <- lapply(X = params, ..., function(i, ..., evalAucPRLoss = NA) {
        OUTRIDER:::evalAutoCorrection(
            ods,
            encoding_dim = i,
            BPPARAM = BPPARAM,
            ...
        )
    })
    metadata(ods)[["encDimTable"]] <- data.table(
        encodingDimension = params,
        evaluationLoss = unlist(eval),
        evalMethod = "aucPR"
    )
    metadata(ods)[["optimalEncDim"]] <- NULL
    metadata(ods)[["optimalEncDim"]] <- OUTRIDER:::getBestQ(ods)
    counts(ods) <- assay(ods, "trueCounts")
    OUTRIDER:::validateOutriderDataSet(ods)
    return(ods)
}

q_params <- round(seq(2, round(ncol(ods)/3), length.out=10))
message('encoding dimension values to search over: ', paste(q_params, collapse = ' '))

# outer_threads <- min(length(q_params), total_threads)
# inner_threads <- max(1, floor(total_threads / outer_threads))

register(MulticoreParam(total_threads))
ods <- findEncodingDim_opt(ods, params = q_params)

saveRDS(ods, snakemake@output$ods)
