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

chooseNumber <- function(observed, truth, max.rank=50)
# Assessing each strategy to choose the number of PCs.
{ 
    center <- rowMeans(observed)
    SVD <- svd(t(observed - center), nu=max.rank, nv=max.rank)
    prog.var <- SVD$d^2 / (ncol(observed) - 1) 
    total.var <- sum(prog.var)

    # Using our denoising approach.
    tech.comp <- apply(observed - truth, 1, var)
    tech.var <- sum(tech.comp)
    denoised <- scran:::.get_npcs_to_keep(prog.var, tech.var)
    
    # Applying parallel analysis.
    parallel <- parallelPCA(observed, BPPARAM=MulticoreParam(3), value="n", approximate=TRUE, min.rank=1, max.rank=max.rank)

    # Using the Marchenko-Pastur limit (following Mathematica's advice).
    library(RMTstat)
    ndf <- nrow(observed) - 1
    limit <- qmp(1, ndf=ndf, pdim=ncol(observed)) * mean(tech.comp)
    marchenko <- sum(SVD$d^2/ndf > limit)

    # Applying the Gavish-Donoho method.
    m <- min(dim(observed))
    n <- max(dim(observed))
    beta <- m/n
    lambda <- sqrt( 2 * (beta + 1) + (8 * beta) / ( beta + 1 + sqrt(beta^2 + 14 * beta + 1) ) )
    gv <- sum(SVD$d > lambda * sqrt(n) * mean(tech.comp))

    # Determining the MSE at each step.
    mse <- computeMSE(SVD, center, truth, ncomponents=max.rank)
    optimal <- which.min(mse)

    return(list(MSE=mse, retained=data.frame(denoised=denoised, parallel=parallel, marchenko=marchenko, gavish=gv, optimal=optimal)))
}

runSimulation <- function(prefix, truth.FUN, iters=10, observed.FUN=NULL)
# A convenience function to run simulations based on a function that generates a matrix of true signal.
{
    fname <- paste0(prefix, ".txt")
    statistics <- list()
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

                is.first <- counter==1L
                write.table(data.frame(Ncells=ncells, Ngenes=ngenes, Prop.DE=affected, cur.retained/iters), 
                    file=fname, append=!is.first, col.names=is.first, row.names=FALSE, quote=FALSE, sep="\t")

                statistics[[counter]] <- cur.stats/iters
                counter <- counter + 1L
            }
        }
    }

    saveRDS(file=paste0(fname, ".rds"), statistics)
    return(NULL)
}

addNoise <- function(variance) 
# Adding noise with different strategies, with different variation of technical noise across genes.
{
    if (variance=="none") {
        function(truth) truth + rnorm(length(truth))
    } else if (variance=="moderate") {
        function(truth) {
            sd <- sqrt(rgamma(nrow(truth), 2, 2)) 
            truth + rnorm(length(truth), sd=sd)
        }
    } else if (variance=="high") {
        function(truth) {
            sd <- sqrt(runif(nrow(truth), 0, 6))
            truth + rnorm(length(truth), sd=sd)
        }
    } else {
        stop("unknown variance mode")
    }
}
