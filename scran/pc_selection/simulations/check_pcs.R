# Setting up a function that decides how many PCs to remove.

library(scran)
library(limma)

assessSeperation <- function(n, pcs, ids) 
# Computing the variance explained by every pair of clusters.
{
    Y <- t(pcs[,seq_len(n),drop=FALSE])
    all.pairs <- combn(length(unique(ids)),2)
   
    collected <- numeric(ncol(all.pairs))
    for (i in seq_along(collected)) { 
        curY <- Y[,ids %in% all.pairs[,i],drop=FALSE]
        rss0 <- sum((curY - rowMeans(curY))^2)
        Y1 <- Y[,ids==all.pairs[1,i],drop=FALSE]
        rss1 <- sum((Y1 - rowMeans(Y1))^2)
        Y2 <- Y[,ids==all.pairs[2,i],drop=FALSE]
        rss2 <- sum((Y2 - rowMeans(Y2))^2)
        collected[i] <- 1 - (rss1 + rss2)/rss0
    }

    return(c(Worst=min(collected), Median=median(collected)))
}

assessPCA <- function(observed, clusters, ids) { 
    x <- prcomp(t(observed)) 
    prog.var <- x$sdev^2 
    total.var <- sum(prog.var)

    truth <- clusters[,ids]
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

    # Assessing each method.
    d.res <- assessSeperation(denoised, x$x, ids)
    p.res <- assessSeperation(parallel, x$x, ids)
    n.res <- assessSeperation(ncol(x$x), x$x, ids)

    return(list(N.denoised=denoised, Worst.denoised=d.res["Worst"], Median.denoised=d.res["Median"],
                N.parallel=parallel, Worst.parallel=p.res["Worst"], Median.parallel=p.res["Median"],
                N.none=ncol(x$x), Worst.none=n.res["Worst"], Median.none=n.res["Median"]))
}

out.tab <- "results_pca.txt"
new.file <- TRUE
set.seed(191919191)

for (ncells in c(200, 1000)) {
    for (ngenes in c(1000, 5000)) {
        for (affected in c(0.2, 0.5, 1)) { 
            for (npops in c(5, 10, 20)) { 
                for (sd in c(0.2, 0.5, 1)) { 
                    
                    collected <- vector("list", 10)
                    for (it in 1:10) {
                        pops <- matrix(rnorm(npops * ngenes, sd=sd), ncol=npops)
                        pops[-seq_len(affected*ngenes),] <- 0

                        in.pop <- sample(npops, ncells, replace=TRUE)
                        true.means <- pops[,in.pop,drop=FALSE]
                        y <- matrix(rnorm(ngenes*ncells, mean=true.means), ncol=ncells)

                        collected[[it]] <- assessPCA(y, pops, in.pop)
                    }

                    collected <- lapply(do.call(mapply, c(collected, list(FUN=c, SIMPLIFY=FALSE))), mean)
                    write.table(file=out.tab, data.frame(Ncells=ncells, Ngenes=ngenes, Prop.DE=affected, 
                                Npops=npops, Bio.SD=sd, collected), sep="\t", quote=FALSE, 
                                row.names=FALSE, col.names=new.file, append=!new.file)
                    new.file <- FALSE
                }
            }
        }
    }
}
