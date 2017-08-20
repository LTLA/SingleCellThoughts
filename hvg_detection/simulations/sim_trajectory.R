# This script simulates trajectory HVGs and compares the performance of the
# two HVG detection methods in scran for detecting trajectories.

library(scran)
source("functions.R")

set.seed(200)
output.tab <- "results_trajectory.txt"
new.file <- TRUE
ngenes <- 1000
ncells <- 500
nspikes <- 50

for (grad in c(-4, -3, -2)) {
    scaling <- matrix(seq_len(ncells) * 10^grad, ngenes, ncells, byrow=TRUE)
    chosen.genes <- 1:50
    scaling[- chosen.genes,] <- 1

    collected <- vector("list", 20)
    for (it in seq_len(20)) {
        sim.genes <- sampleCounts(ngenes=ngenes, nsamples=ncells, scaling=scaling)
        sim.spikes <- sampleCounts(ngenes=nspikes, nsamples=ncells, scaling=1)
        counts <- rbind(sim.genes, sim.spikes)
        is.spike <- seq_len(nspikes) + ngenes
        my.env <- detectHVGs(counts, is.spike)
        collected[[it]] <- my.env$output
    }

    write.table(data.frame(Grad=grad, rbind(colMeans(do.call(rbind, collected)))),
                file=output.tab, append=!new.file, col.names=new.file, row.names=FALSE, sep="\t", quote=FALSE)
    new.file <- FALSE
}
