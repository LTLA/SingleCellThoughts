# This script generates gaussian cluster simulations.
# It expects:
# - NPOPS: the number of populations.
# - SD: the standard deviation of the population means.

if (!exists(NPOPS)) NPOPS <- 5
if (!exists(SD)) SD <- 1

FUN <- function(ngenes, ncells) {
    in.pop <- sample(NPOPS, ncells, replace=TRUE)
    pop.mean <- matrix(rnorm(NPOPS * ngenes, sd=SD), ncol=NPOPS)
    pop.mean[,in.pop,drop=FALSE]
}

source("functions.R")
dir.create("results", showWarnings=FALSE)
runSimulation(paste0("results/gaussclust-%s-%s.txt", NPOPS, SD), FUN)
