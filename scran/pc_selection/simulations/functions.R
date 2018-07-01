# Setting up a function that decides how many PCs to remove.

library(scran)

computeAllMetrics <- function(truth, iterations=200, upper=50) {
    fitted <- as.list(numeric(upper))
    mse <- numeric(upper)

    for (it in seq_len(iterations)) {
        observed <- truth + rnorm(length(truth))
        center <- rowMeans(observed)
        svd.out <- svd(t(observed - center), nu=upper, nv=upper)
        running.recon <- 0

        for (N in seq_len(upper)) {
            current <- outer(svd.out$v[,N] * svd.out$d[N], svd.out$u[,N])
            running.recon <- running.recon + current
            full.recon <- running.recon + center 
        
            mse[N] <- mse[N] + mean((full.recon - truth)^2)
            fitted[[N]] <- fitted[[N]] + full.recon
        }
    }

    mse <- mse/iterations
    bias <- numeric(upper)
    for (N in seq_len(upper)) {
        fitted[[N]] <- fitted[[N]]/iterations
        bias[N] <- mean((fitted[[N]] - truth)^2)
    }
    
    return(data.frame(Bias=bias, Variance=mse-bias, MSE=mse))
}

chooseNumber <- function(observed, truth)
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
    upper <- max(1L, sum(attr(parallel, "percentVar") > top.var))

    return(data.frame(denoised=denoised, parallel=parallel, upper=upper))
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
                cur.stats <- cur.retained <- NULL 

                for (it in seq_len(iters)) {
                    truth <- FUN(ngenes*affected, ncells)
                    truth <- rbind(truth, matrix(0, ncol=ncells, nrow=(1-affected)*ngenes))
                    metrics <- computeAllMetrics(truth)

                    y <- truth + rnorm(length(truth))
                    retained <- chooseNumber(y, truth)

                    if (it==1L) {
                        cur.stats <- metrics
                        cur.retained <- retained
                    } else {
                        cur.stats <- cur.stats + metrics
                        cur.retained <- cur.retained + retained
                    }
                }

                cur.retained <- cur.retained/iters
                cur.stats <- cur.stats/iters

                scenarios[[counter]] <- data.frame(Ncells=ncells, Ngenes=ngenes, Prop.DE=affected)
                statistics[[counter]] <- cur.stats
                numbers[[counter]] <- cur.retained
                counter <- counter + 1L
            }
        }
    }

    scenarios <- do.call(rbind, scenarios)
    numbers <- do.call(rbind, numbers)
    saveRDS(file=fname, list(scenarios=scenarios, statistics=statistics, retained=numbers))
    return(NULL)
}
