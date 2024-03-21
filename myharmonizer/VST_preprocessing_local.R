
suppressPackageStartupMessages(library("DESeq2"))

## Modified functions from DESeq2 #####################################################################################################
# Original functions available under estimateSizeFactorsForMatrix_DESeq2.R or at https://code.bioconductor.org/browse/DESeq2/blob/RELEASE_3_12/R/core.R
# Modified from DESeq2 Release_3_12 (19 Feb 2021)
# Citation: Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. \url{https://doi.org/10.1186/s13059-014-0550-8}  
# Last modified on 24 November 2023.

estimateSizeFactorsForMatrix.modified <- function (counts, locfunc = stats::median, geoMeans, controlGenes,
                                                   type = c("ratio", "poscounts"))
{
  type <- match.arg(type, c("ratio", "poscounts"))
  if (missing(geoMeans)) {
    incomingGeoMeans <- FALSE
    if (type == "ratio") {
      loggeomeans <- rowMeans(log(counts))
    }
    else if (type == "poscounts") {
      lc <- log(counts)
      lc[!is.finite(lc)] <- 0
      loggeomeans <- rowMeans(lc)
      allZero <- rowSums(counts) == 0
      loggeomeans[allZero] <- -Inf
    }
  }
  else {
    incomingGeoMeans <- TRUE
    if (length(geoMeans) != nrow(counts)) {
      stop("geoMeans should be as long as the number of rows of counts")
    }
    loggeomeans <- log(geoMeans)
  }
  if (all(is.infinite(loggeomeans))) {
    stop("every gene contains at least one zero, cannot compute log geometric means")
  }
  sf <- if (missing(controlGenes)) {
    apply(counts, 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) &
                                              cnts > 0]))
    })
  }
  else {
    if (!(is.numeric(controlGenes) | is.logical(controlGenes))) {
      stop("controlGenes should be either a numeric or logical vector")
    }
    loggeomeansSub <- loggeomeans[controlGenes]
    apply(counts[controlGenes, , drop = FALSE], 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeansSub)[is.finite(loggeomeansSub) &
                                                 cnts > 0]))
    })
  }

  ## Modified from DESeq2 Release_3_12 (19 Feb 2021).
  ## Last modified on 24 November 2023. Modification was to comment out two lines below so that size
  ## factors are not stabilized in the frozen representation (as this unfreezes the transformations).
  # if (incomingGeoMeans) {
  #   sf <- sf/exp(mean(log(sf)))
  # }
  sf
}

## End Modified functions from DESeq2 #################################################################################################

vst_preprocessing <- function(data_r, geneCorr = c('none', 'rpk'), dispersionList=NULL, output=getwd()){

  geneCorr <- match.arg(geneCorr)

  if(is.null(dispersionList)){
    print('myHarmonizer object does not have the frozen dispersion function.')

  } else {
    switch(geneCorr,
           'none' = {
             data <- data_r
             prefix = 'VST'
           },
           'rpk' = {

             #import genelength
             # genelength = read.csv(geneLengthFile)
             genelength = dispersionList$genelength

             #preprocess with gene length
             data = data_r[,colnames(data_r) %in% genelength$gene]
             genelength_g = genelength$med[match(colnames(data), genelength$gene)]

             if(ncol(data) == length(genelength_g)){
               data = ceiling(sweep(data, 2, as.array(genelength_g), `/`)*10000)
             } else {stop("Mismatch between number of genes in sample and number of genes in genelength variable.")}

             prefix = 'GeVST'
           })

    # #track and remove genes with zero expression
    # zeroLoc = (apply(data, 2, sum) == 0)
    # dataz = data[,!zeroLoc]

    colData = data.frame('sample'=rownames(data))
    dds = DESeqDataSetFromMatrix(countData = t(data),
                                 colData = colData,
                                 design = ~1)


    sizeFactors(dds) <- estimateSizeFactorsForMatrix.modified(counts(dds), geoMeans = dispersionList$gm)
    dispersionFunction(dds) <- dispersionList$dispersion

    ####WARNING!!! Size factors will no longer multiply to 1 for validation, test datasets

    vst = varianceStabilizingTransformation(dds, blind=FALSE) #False in order to use dispersion estimates from DESeq

    vst_mz = t(assay(vst))

    return(vst_mz)

  }
}
