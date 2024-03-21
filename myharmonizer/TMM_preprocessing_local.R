suppressPackageStartupMessages(library("edgeR"))

## Modified functions from EdgeR #####################################################################################################
# Original functions available under calcNormFactors_edgeR.R or at https://code.bioconductor.org/browse/edgeR/blob/RELEASE_3_16/R/calcNormFactors.R
# Modified from edgeR Release_3_16 (modified 2 June 2020)
# Citation:   Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential
              # expression analysis of digital gene expression data. Bioinformatics 26, 139-140

              # McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression analysis of multifactor RNA-Seq
              # experiments with respect to biological variation. Nucleic Acids Research 40, 4288-4297

              # Chen Y, Lun ATL, Smyth GK (2016). From reads to genes to pathways: differential expression
              # analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
              # F1000Research 5, 1438

              # Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK (2024). edgeR 4.0: powerful differential analysis
              # of sequencing data with expanded functionality and improved support for small counts and larger
              # datasets. bioRxiv doi: 10.1101/2024.01.21.576131
# Last modified on 1 November 2023.

calcNormFactors.modified <- function(object, lib.size=NULL, method=c("TMM","TMMwsp","RLE","upperquartile","none"), refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, factor_scaling=NULL, ...)
  #	Scale normalization of RNA-Seq data, for count matrices
  #	Mark Robinson, Gordon Smyth and edgeR team
  #	Created 22 October 2009. Last modified 2 June 2020.
{
  #	Check object
  x <- as.matrix(object)
  if(any(is.na(x))) stop("NA counts not permitted")
  nsamples <- ncol(x)

  #	Check lib.size
  if(is.null(lib.size)) {
    lib.size <- colSums(x)
  } else {
    if(anyNA(lib.size)) stop("NA lib.sizes not permitted")
    if(length(lib.size) != nsamples) {
      if(length(lib.size) > 1L) warning("calcNormFactors: length(lib.size) doesn't match number of samples",call.=FALSE)
      lib.size <- rep_len(lib.size,nsamples)
    }
  }

  #	Check method
  #	Backward compatability with previous name
  if(length(method)==1L && method=="TMMwzp") {
    method <- "TMMwsp"
    message("TMMwzp has been renamed to TMMwsp")
  }
  method <- match.arg(method)

  #	Remove all zero rows
  allzero <- .rowSums(x>0, nrow(x), nsamples) == 0L
  if(any(allzero)) x <- x[!allzero,,drop=FALSE]

  #	Degenerate cases
  if(nrow(x)==0 || nsamples==1) method="none"

  #	Calculate factors
  f <- switch(method,
              TMM = {
                if( is.null(refColumn) ) {
                  f75 <- suppressWarnings(.calcFactorQuantile(data=x, lib.size=lib.size, p=0.75))
                  if(median(f75) < 1e-20) {
                    refColumn <- which.max(colSums(sqrt(x)))
                  } else {
                    refColumn <- which.min(abs(f75-mean(f75)))
                  }
                }
                f <- rep_len(NA_real_,nsamples)
                for(i in 1:nsamples)
                  f[i] <- .calcFactorTMM(obs=x[,i],ref=x[,refColumn], libsize.obs=lib.size[i], libsize.ref=lib.size[refColumn], logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
                f
              },
              TMMwsp = {
                if( is.null(refColumn) ) refColumn <- which.max(colSums(sqrt(x)))
                f <- rep_len(NA_real_,nsamples)
                for(i in 1:nsamples)
                  f[i] <- .calcFactorTMMwsp(obs=x[,i],ref=x[,refColumn], libsize.obs=lib.size[i], libsize.ref=lib.size[refColumn], logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
                f
              },
              RLE = .calcFactorRLE(x)/lib.size,
              upperquartile = .calcFactorQuantile(x,lib.size,p=p),
              none = rep_len(1,nsamples)
  )

  #	Factors should multiple to one
  ## Modified from edgeR Release_3_16 modified 2 June 2020.
  ## Last modified on 1 November 2023.
  if(is.null(factor_scaling)) factor_scaling = exp(mean(log(f)))
  f <- f/factor_scaling

  #	Output
  ## Modified from edgeR Release_3_16 modified 2 June 2020.
  ## Last modified on 1 November 2023.
  names(f) <- colnames(x)
  list(f = f,
       ref = refColumn,
       factor_scaling = factor_scaling)
}

.calcFactorRLE <- function(data)
  #	Scale factors as in Anders et al (2010)
  #	Mark Robinson
  #	Created 16 Aug 2010
{
  gm <- exp(rowMeans(log(data)))
  apply(data, 2, function(u) median((u/gm)[gm > 0]))
}

.calcFactorQuantile <- function (data, lib.size, p=0.75)
  #	Generalized version of upper-quartile normalization
  #	Mark Robinson and Gordon Smyth
  #	Created 16 Aug 2010. Last modified 12 Sep 2020.
{
  f <- rep_len(1,ncol(data))
  for (j in seq_len(ncol(data))) f[j] <- quantile(data[,j], probs=p)
  if(min(f)==0) warning("One or more quantiles are zero")
  f / lib.size
}

.calcFactorTMM <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
  #	TMM between two libraries
  #	Mark Robinson
{
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)

  if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs
  if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref

  logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance

  #	remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]

  if(max(abs(logR)) < 1e-6) return(1)

  #	taken from the original mean() function
  n <- length(logR)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS

  #	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  #	a fix from leonardo ivan almonacid cardenas, since rank() can return
  #	non-integer values when there are a lot of ties
  keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

  if(doWeighting)
    f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
  else
    f <- mean(logR[keep], na.rm=TRUE)

  #	Results will be missing if the two libraries share no features with positive counts
  #	In this case, return unity
  if(is.na(f)) f <- 0
  2^f
}

.calcFactorTMMwsp <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
  #	TMM with pairing of singleton positive counts between the obs and ref libraries
  #	Gordon Smyth
  #	Created 19 Sep 2018. Last modified 9 Jun 2020.
{
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)

  #	epsilon serves as floating-point zero
  eps <- 1e-14

  #	Identify zero counts
  pos.obs <- (obs > eps)
  pos.ref <- (ref > eps)
  npos <- 2L * pos.obs + pos.ref

  #	Remove double zeros and NAs
  i <- which(npos==0L | is.na(npos))
  if(length(i)) {
    obs <- obs[-i]
    ref <- ref[-i]
    npos <- npos[-i]
  }

  #	Check library sizes
  if(is.null(libsize.obs)) libsize.obs <- sum(obs)
  if(is.null(libsize.ref)) libsize.ref <- sum(ref)

  #	Pair up as many singleton positives as possible
  #	The unpaired singleton positives are discarded so that no zeros remain
  zero.obs <- (npos == 1L)
  zero.ref <- (npos == 2L)
  k <- (zero.obs | zero.ref)
  n.eligible.singles <- min( sum(zero.obs), sum(zero.ref))
  if(n.eligible.singles > 0L) {
    refk <- sort(ref[k],decreasing=TRUE)[1:n.eligible.singles]
    obsk <- sort(obs[k],decreasing=TRUE)[1:n.eligible.singles]
    obs <- c(obs[!k],obsk)
    ref <- c(ref[!k],refk)
  } else {
    obs <- obs[!k]
    ref <- ref[!k]
  }

  #	Any left?
  n <- length(obs)
  if(n==0L) return(1)

  #	Compute M and A values
  obs.p <- obs / libsize.obs
  ref.p <- ref / libsize.ref
  M <- log2( obs.p / ref.p )
  A <- 0.5 * log2( obs.p * ref.p )

  #	If M all zero, return 1
  if(max(abs(M)) < 1e-6) return(1)

  #	M order, breaking ties by shrunk M
  obs.p.shrunk <- (obs+0.5) / (libsize.obs+0.5)
  ref.p.shrunk <- (ref+0.5) / (libsize.ref+0.5)
  M.shrunk <- log2( obs.p.shrunk / ref.p.shrunk )
  o.M <- order(M, M.shrunk)

  #	A order
  o.A <- order(A)

  #	Trim
  loM <- as.integer(n * logratioTrim) + 1L
  hiM <- n + 1L - loM
  keep.M <- rep_len(FALSE,n)
  keep.M[o.M[loM:hiM]] <- TRUE
  loA <- as.integer(n * sumTrim) + 1L
  hiA <- n + 1L - loA
  keep.A <- rep_len(FALSE,n)
  keep.A[o.A[loA:hiA]] <- TRUE
  keep <- keep.M & keep.A
  M <- M[keep]

  #	Average the M values
  if(doWeighting) {
    obs.p <- obs.p[keep]
    ref.p <- ref.p[keep]
    v <- (1-obs.p)/obs.p/libsize.obs + (1-ref.p)/ref.p/libsize.ref
    w <- (1+1e-6) / (v+1e-6)
    TMM <- sum(w*M) / sum(w)
  } else {
    TMM <- mean(M)
  }

  2^TMM
 }

## End Modified functions from EdgeR #################################################################################################

tmm_preprocessing <- function(data, referenceSample=NULL){

  data = t(data)

  #determine if referenceSample is null (training dataset)
  if(is.null(referenceSample)){

    print('myHarmonizer object does not have the frozen reference sample.')

  } else {
    refs = referenceSample
    data_f = cbind(refs$ref, data)
    data_list_tmm = calcNormFactors.modified(data_f, refColumn = 1, factor_scaling = refs$factor_scaling)

    if(ncol(data_f)!=length(data_list_tmm$f)){stop("Mismatch between number of samples and number of scaling factors.")
    } else { #apply scaling factor
      data_scaled_plus = sweep(data_f, 2, as.array(data_list_tmm$f), `*`)
      #log cpm transformation (EdgeR cpm preprocessing)
      data_cpm_plus = log(sapply(seq(ncol(data_f)), function(c){
        (data_f[,c]/sum(data_scaled_plus[,c]))*1e6}))
      colnames(data_cpm_plus) <- colnames(data_scaled_plus)

      return(t(data_cpm_plus[,-1]))
    }
  }
}

g_tmm_preprocessing <- function(data, gene_length=NULL, referenceSample=NULL){

  data = t(data)

  if(is.null(referenceSample)){
    print('myHarmonizer object does not have the frozen reference sample.')

  } else {
    refs = referenceSample
    genelength = refs$genelength
  }

  #preprocess with gene length
  data_g = data[rownames(data) %in% genelength$gene,]
  genelength_g = genelength$med[match(rownames(data_g), genelength$gene)]

  if(nrow(data_g) == length(genelength_g)){
    data_g = sweep(data_g, 1, as.array(genelength_g), `/`)
  } else {stop("Mismatch between number of genes in sample and number of genes in genelength variable.")}

      data_fg = cbind(refs$ref_g, data_g)

    data_list_tmm_g = calcNormFactors.modified(data_fg, refColumn = 1, factor_scaling = refs$factor_scaling_g)

    if(ncol(data_fg)!=length(data_list_tmm_g$f)){stop("Mismatch between number of samples and number of scaling factors.")
    } else {
      data_g_scaled_plus = sweep(data_fg, 2, as.array(data_list_tmm_g$f), `*`)
      data_g_cpm_plus = log(sapply(seq(ncol(data_fg)), function(c){
        (data_fg[,c]/sum(data_g_scaled_plus[,c]))*1e6}))
      colnames(data_g_cpm_plus) <- colnames(data_g_scaled_plus)

      return(t(data_g_cpm_plus[,-1]))
    }
}

