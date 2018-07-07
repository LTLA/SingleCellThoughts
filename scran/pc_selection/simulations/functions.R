# Setting up a function that decides how many PCs to remove.

library(scran)

computeMSE <- function(svd.out, center, truth, ncomponents=100) {
    running <- 0
    collected <- numeric(ncomponents)
    for (i in seq_along(collected)) { 
        running <- running + svd.out$u[,i,drop=FALSE] %*% (svd.out$d[i] * t(svd.out$v[,i,drop=FALSE]))
        collected[i] <- mean((t(running) + center - truth)^2)
    }
    return(collected)
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
    tech.comp <- apply(observed - truth, 1, var)
    tech.var <- sum(tech.comp)
    denoised <- scran:::.get_npcs_to_keep(prog.var, tech.var)
    
    # Applying parallel analysis.
    parallel <- parallelPCA(observed, BPPARAM=MulticoreParam(3), value="n", approximate=TRUE, min.rank=1)

    # A quick-and-dirty threshold on the upper bound on random eigenvalues.
    top.var <- colMeans(attr(parallel, "permuted.percentVar"))[1]
    upper <- sum(attr(parallel, "percentVar") > top.var)

    # Applying the Gavish-Donoho method.
    m <- min(dim(observed))
    n <- max(dim(observed))
    beta <- m/n
    lambda <- sqrt( 2 * (beta + 1) + (8 * beta) / ( beta + 1 + sqrt(beta^2 + 14 * beta + 1) ) )
    gv <- sum(SVD$d > lambda * sqrt(n) * mean(tech.comp))

    # Determining the MSE at each step.
    mse <- computeMSE(SVD, center, truth)
    optimal <- which.min(mse)

    return(list(MSE=mse, retained=data.frame(denoised=denoised, parallel=parallel, upper=upper, gavish=gv, optimal=optimal)))
}

runSimulation <- function(fname, truth.FUN, iters=10, observed.FUN=NULL)
# A convenience function to run simulations based on a function that generates a matrix of true signal.
{
    scenarios <- list()
    statistics <- list()
    numbers <- list()
    counter <- 1L

    if (is.null(observed.FUN)) {
        observed.FUN <- function(truth) {
            truth + rnorm(length(truth))
        }
    }

    for (ncells in c(200, 1000)) {
        for (ngenes in c(1000, 5000)) {
            for (affected in c(0.2, 0.5, 1)) { 
                cur.mse <- cur.retained <- NULL 

                for (it in seq_len(iters)) {
                    truth <- truth.FUN(ngenes*affected, ncells)
                    truth <- rbind(truth, matrix(0, ncol=ncells, nrow=(1-affected)*ngenes))
                    y <- observed.FUN(truth)
                    out <- chooseNumber(y, truth)

                    if (it==1L) {
                        cur.stats <- out$MSE
                        cur.retained <- out$retained
                    } else {
                        cur.stats <- cur.stats + out$MSE
                        cur.retained <- cur.retained + out$retained
                    }
                }

                cur.stats <- cur.stats/iters
                cur.retained <- cur.retained/iters

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
