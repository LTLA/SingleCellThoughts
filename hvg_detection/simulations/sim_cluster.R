# This script simulates cluster HVGs and compares the performance of the
# two HVG detection methods in scran for detecting clusters.

library(scran)
source("functions.R")

set.seed(200)
output.tab <- "results_cluster.txt"
new.file <- TRUE
ngenes <- 1000
ncells <- 500
nspikes <- 50

for (subpop in c(10, 20, 50, 100, 200)) { 
    chosen.cells <- seq_len(subpop)

    for (fold in c(2, 5, 10)) {
        scaling <- matrix(1, ngenes, ncells)
        chosen.genes <- 1:50
        scaling[chosen.genes, chosen.cells] <- fold

        collected <- vector("list", 20)
        for (it in seq_len(20)) {
            sim.genes <- sampleCounts(ngenes=ngenes, nsamples=ncells, scaling=scaling)
            sim.spikes <- sampleCounts(ngenes=nspikes, nsamples=ncells, scaling=1)
            counts <- rbind(sim.genes, sim.spikes)
            is.spike <- seq_len(nspikes) + ngenes
            my.env <- detectHVGs(counts, is.spike)
            collected[[it]] <- my.env$output
        }

        write.table(data.frame(Ncells=subpop, FC=fold, rbind(colMeans(do.call(rbind, collected)))),
                    file=output.tab, append=!new.file, col.names=new.file, row.names=FALSE, sep="\t", quote=FALSE)
        new.file <- FALSE
    }
}

