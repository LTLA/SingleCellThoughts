# Setting up a function that decides how many PCs to remove.

library(scran)
library(limma)

assessAccuracy <- function(n, svd.out, center, truth)
# Computing the MSE of the reconstruction with a given number of PCs. 
{
    out <- svd.out$u[,1:n] %*% (svd.out$d[1:n] * t(svd.out$v[,1:n])) 
    out <- t(out) + center
    return(sum((out-truth)^2))
}

assessPCA <- function(observed, truth) 
# Assessing each strategy to choose the number of PCs.
{ 
    center <- rowMeans(observed)
    x <- svd(t(observed - center))
    prog.var <- x$d^2 / (ncol(observed) - 1) 
    total.var <- sum(prog.var)

    tech.var <- sum(apply(observed - truth, 1, var))
    npcs <- length(prog.var)
    flipped.prog.var <- rev(prog.var)

## Too unstable to fluctuations!
#    estimated.contrib <- cumsum(flipped.prog.var) + flipped.prog.var * (npcs:1 - 1L)
#    estimated.contrib <- rev(estimated.contrib)
#    retained <- min(which(tech.var > estimated.contrib))

    # Using our denoising approach, or parallel analysis.
    denoised <- scran:::.get_npcs_to_keep(prog.var, tech.var)
    approximate <- ncol(observed)>=500
    parallel <- parallelPCA(observed, BPPARAM=MulticoreParam(3), value="n", approximate=approximate, min.rank=1)

    # A quick-and-dirty Marchenko-Pastur-like threshold.
    marpast <- sum(attr(parallel, "percentVar") > colMeans(attr(parallel, "permuted.percentVar"))[1])

    # Assessing each method.
    d.res <- assessAccuracy(denoised, x, center, truth)
    p.res <- assessAccuracy(parallel, x, center, truth)
    m.res <- assessAccuracy(marpast, x, center, truth)
    n.res <- assessAccuracy(ncol(x$u), x, center, truth)

    return(list(N.denoised=denoised, MSE.denoised=d.res,
                N.parallel=parallel, MSE.parallel=p.res,
                N.marpast=marpast, MSE.marpast=m.res,
                N.none=ncol(x$x), MSE.none=n.res))
}

runSimulation <- function(fname, FUN) 
# A convenience function to run simulations based on a function that 
# generates a matrix of true signal.
{
    new.file <- TRUE
    for (ncells in c(200, 1000)) {
        for (ngenes in c(1000, 5000)) {
            for (affected in c(0.2, 0.5, 1)) { 

                collected <- vector("list", 10)
                for (it in 1:10) {
                    truth <- FUN(ngenes, ncells)
                    truth[-seq_len(affected*ngenes),] <- 0
                    y <- matrix(rnorm(ngenes*ncells, mean=truth), ncol=ncells)
                    collected[[it]] <- assessPCA(y, truth)
                }

                collected <- lapply(do.call(mapply, c(collected, list(FUN=c, SIMPLIFY=FALSE))), mean)
                write.table(file=fname, data.frame(Ncells=ncells, Ngenes=ngenes, Prop.DE=affected, collected), 
                    sep="\t", quote=FALSE, row.names=FALSE, col.names=new.file, append=!new.file)
                new.file <- FALSE
            }
        }
    }
    return(NULL)
}
