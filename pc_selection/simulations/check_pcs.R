# Setting up a function that decides how many PCs to remove.

assessPCA <- function(observed, truth) { 
    x <- prcomp(t(observed)) 
    prog.var <- x$sdev^2 
    total.var <- sum(prog.var)

## Sanity checks.    
#    stopifnot(isTRUE(all.equal(total.var, sum(apply(observed, 1, var)))))
#    stopifnot(isTRUE(all.equal(total.var, sum(apply(x$x, 2, var)))))
    
    tech.var <- sum(apply(observed - truth, 1, var))
    npcs <- length(prog.var)
    flipped.prog.var <- rev(prog.var)

## Too unstable to fluctuations!
#    estimated.contrib <- cumsum(flipped.prog.var) + flipped.prog.var * (npcs:1 - 1L)
#    estimated.contrib <- rev(estimated.contrib)
#    retained <- min(which(tech.var > estimated.contrib))
    retained <- npcs - min(which(cumsum(prog.var) > tech.var)) + 1

    # These aren't quite right, they only approximate the two components.
    ps.bio.pos <- t(truth - x$center) %*% x$rotation 
    ps.bio.var <- apply(ps.bio.pos, 2, var)

    scatter <- observed - truth
    ps.tech.pos <- t(scatter) %*% x$rotation 
    ps.tech.var <- apply(ps.tech.pos, 2, var)

    return(list(bio=ps.bio.var, tech=ps.tech.var, total=prog.var, retained=retained))
}

out.tab <- "results_pca.txt"
new.file <- TRUE

for (ncells in c(200, 500, 1000)) {
    for (ngenes in c(1000, 2000, 5000)) {
        for (affected in c(0.2, 0.5, 1)) { 
            for (npops in c(5, 10, 20)) { 
                for (sd in c(0.2, 0.5, 1)) { 
                    
                    nkept <- tech.prop <- bio.prop <- old.kept <- vector("list", 10)
                    for (it in 1:10) {
                        pops <- matrix(rnorm(npops * ngenes, sd=sd), ncol=npops)
                        in.pop <- sample(npops, ncells, replace=TRUE)
                        true.means <- pops[,in.pop,drop=FALSE]
                        true.means[-seq_len(affected*ngenes),] <- 0
                        y <- matrix(rnorm(ngenes*ncells, mean=true.means), ncol=ncells)
                       
                        observed<-y
                        truth <- true.means
                        out <- assessPCA(y, true.means)
                        bio.prop[[it]] <- sum(out$bio[seq_len(out$retained)])/sum(out$bio)
                        tech.prop[[it]] <- sum(out$tech[seq_len(out$retained)])/sum(out$tech)
                        nkept[[it]] <- out$retained
                    }

                    nkept <- mean(unlist(nkept))
                    tech.prop <- mean(unlist(tech.prop))
                    bio.prop <- mean(unlist(bio.prop))
                    write.table(file=out.tab, data.frame(Ncells=ncells, Ngenes=ngenes, Prop.DE=affected,
                                                         Npops=npops, Bio.SD=sd, NPC=nkept, 
                                                         Bio.Prop=bio.prop, Tech.Prop=tech.prop),
                                sep="\t", quote=FALSE, row.names=FALSE, col.names=new.file, append=!new.file)
                    new.file <- FALSE
                }
            }
        }
    }
}



