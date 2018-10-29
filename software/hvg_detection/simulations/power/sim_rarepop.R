# This simulates HVGs associated with rare populations.

library(scran)
source("../functions.R")

set.seed(200)
detect.tab <- "detect_rarepop.txt"
above.tab <- "above_rarepop.txt"
new.file <- TRUE
ngenes <- 1000
ncells <- 500
nspikes <- 50

subpop <- 5
for (scenario in c("up", "down")) { 
    for (fold in c(0.5, 1, 2, 5)) { 
        scaling <- matrix(1, ngenes, ncells)
        chosen.genes <- 1:50

        # Either genes are upregulated or downregulated in the rare population.
        if (scenario=="up") {
            scaling[chosen.genes, seq_len(subpop)] <- fold
            scaling[chosen.genes, -seq_len(subpop)] <- 0
        } else {
            scaling[chosen.genes, seq_len(subpop)] <- 0
            scaling[chosen.genes, -seq_len(subpop)] <- fold
        }

        collected.d <- collected.a <- vector("list", 20)
        for (it in seq_len(20)) {
            sim.genes <- sampleCounts(ngenes=ngenes, nsamples=ncells, scaling=scaling)
            sim.spikes <- sampleCounts(ngenes=nspikes, nsamples=ncells, scaling=1)
            counts <- rbind(sim.genes, sim.spikes)
            is.spike <- seq_len(nspikes) + ngenes

            my.env <- detectHVGs(counts, is.spike, chosen.genes)
            collected.d[[it]] <- my.env$output
            collected.a[[it]] <- my.env$above
        }

        write.table(data.frame(Direction=scenario, FC=fold, rbind(colMeans(do.call(rbind, collected.d)))),
                    file=detect.tab, append=!new.file, col.names=new.file, row.names=FALSE, sep="\t", quote=FALSE)
        write.table(data.frame(Direction=scenario, FC=fold, rbind(colMeans(do.call(rbind, collected.a)))),
                    file=above.tab, append=!new.file, col.names=new.file, row.names=FALSE, sep="\t", quote=FALSE)
        new.file <- FALSE
    }
}




