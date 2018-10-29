computeDisp <- function(means, scale=10, asymptote=0.2) {
    pmin(1000, scale/means + asymptote)
}

sampleCounts <- function(ngenes, nsamples, bounds=c(-2, 10), scaling=1, ...) { 
    means <- 2^runif(ngenes, bounds[1], bounds[2]) * scaling 
    dispersions <- computeDisp(means, ...)
    counts <- matrix(rnbinom(ngenes*nsamples, mu=means, size=1/dispersions), 
                     nrow=ngenes, ncol=nsamples) 
    return(counts)
} 

assessDetection <- function(fdr, to.change, threshold) { 
    is.sig <- which(fdr <= threshold)
    top <- order(fdr)[1:100]
    return(c(True=length(intersect(is.sig, to.change)), 
             Total=length(is.sig),
             Top=length(intersect(top, to.change))))
}

assessAbove <- function(val, trend, to.change) { 
    return(c(Above=sum((val>trend)[to.change], na.rm=TRUE),
             Total=sum(val>trend, na.rm=TRUE)))
}


detectHVGs <- function(counts, is.spike, chosen, threshold=0.05) { 
    # Technical CV2
    out.raw <- technicalCV2(counts, is.spike, sf.cell=1, sf.spike=1, min.bio.disp=0)
    raw.cv2 <- assessDetection(out.raw$FDR, chosen, threshold)
    
    # Variance of log-counts
    lcounts <- log2(counts+1)
    fit <- suppressWarnings(trendVar(lcounts, span=0.2, subset.row=is.spike))
    ref <- decomposeVar(lcounts, fit, test="f")
    raw.log <- assessDetection(ref$FDR, chosen, threshold)
    
    # Improved CV2.
    out.imp <- improvedCV2(counts, is.spike, sf.cell=1, sf.spike=1)
    improved.cv2 <- assessDetection(out.imp$FDR, chosen, threshold)
    
    output <- c(CV2=raw.cv2, Log=raw.log, ImprovedCV2=improved.cv2)
    above <- c(CV2=assessAbove(out.raw$cv2, out.raw$trend, chosen),
               Log=assessAbove(ref$total, ref$tech, chosen),
               ImprovedCV2=assessAbove(out.imp$cv2, out.imp$trend, chosen))
    return(environment())
}
