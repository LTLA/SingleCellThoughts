# Setting up a function that decides how many PCs to remove.

library(scran)
library(limma)

computeAllMetrics <- function(svd.out, center, truth, design, upper=50) {
    bias <- numeric(upper)
    rss <- numeric(upper)
    mse <- numeric(upper)
    running.recon <- 0

    for (N in seq_len(upper)) {
        current <- outer(svd.out$v[,N] * svd.out$d[N], svd.out$u[,N])
        running.recon <- running.recon + current
        full.recon <- running.recon + center  
        fit <- lmFit(full.recon, design)
        bias[N] <- mean((fit$coefficients %*% t(design) - truth)^2)
        rss[N] <- mean(fit$sigma^2)    
        mse[N] <- mean((full.recon - truth)^2)
    }
    
    return(data.frame(Bias=bias, Variance=rss, MSE=mse))
}

assessPCA <- function(observed, truth, design)
# Assessing each strategy to choose the number of PCs.
{ 
    center <- rowMeans(observed)
    SVD <- svd(t(observed - center))
    prog.var <- SVD$d^2 / (ncol(observed) - 1) 
    total.var <- sum(prog.var)

## Too unstable to fluctuations!
#    npcs <- length(prog.var)
#    flipped.prog.var <- rev(prog.var)
#    estimated.contrib <- cumsum(flipped.prog.var) + flipped.prog.var * (npcs:1 - 1L)
#    estimated.contrib <- rev(estimated.contrib)
#    retained <- min(which(tech.var > estimated.contrib))

    # Using our denoising approach.
    tech.var <- sum(apply(observed - truth, 1, var))
    denoised <- scran:::.get_npcs_to_keep(prog.var, tech.var)
    
    # Applying parallel analysis.
    approximate <- ncol(observed)>=500
    parallel <- parallelPCA(observed, BPPARAM=MulticoreParam(3), value="n", approximate=approximate, min.rank=1)

    # A quick-and-dirty threshold on the upper bound on random eigenvalues.
    top.var <- colMeans(attr(parallel, "permuted.percentVar"))[1]
    upper <- sum(attr(parallel, "percentVar") > top.var)

    # Assessing all PC numbers.
    swept <- computeAllMetrics(SVD, center, truth, design)
    optimal <- which.min(swept$MSE)

    return(list(metrics=swept, retained=data.frame(denoised=denoised, parallel=parallel, upper=upper, optimal=optimal)))
}

runSimulation <- function(fname, FUN, iters=10) 
# A convenience function to run simulations based on a function that 
# generates a matrix of true signal.
{
    scenarios <- list()
    statistics <- list()
    numbers <- list()
    counter <- 1L

    for (ncells in c(200, 1000)) {
        for (ngenes in c(1000, 5000)) {
            for (affected in c(0.2, 0.5, 1)) { 
                cur.stats <- cur.n <- NULL 

                for (it in seq_len(iters)) {
                    out <- FUN(ngenes*affected, ncells)
                    truth <- rbind(out$truth, matrix(0, ncol=ncells, nrow=(1-affected)*ngenes))
                    y <- matrix(rnorm(ngenes*ncells, mean=truth), ncol=ncells)
                    assessed <- assessPCA(y, truth, out$design)

                    if (it==1L) {
                        cur.stats <- assessed$metrics
                        cur.n <- assessed$retained
                    } else {
                        cur.stats <- cur.stats + assessed$metrics
                        cur.n <- cur.n + assessed$retained
                    }
                }

                cur.n <- cur.n/iters
                cur.stats <- cur.stats/iters

                scenarios[[counter]] <- data.frame(Ncells=ncells, Ngenes=ngenes, Prop.DE=affected)
                statistics[[counter]] <- cur.stats
                numbers[[counter]] <- cur.n
                counter <- counter + 1L
            }
        }
    }

    saveRDS(file=fname, list(scenarios=scenarios, statistics=statistics, retained=numbers))
    return(NULL)
}
