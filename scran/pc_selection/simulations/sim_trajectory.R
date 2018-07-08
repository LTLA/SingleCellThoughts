# This script generates trajectories between cell populations.
# It expects:
# - NPOPS: the number of populations. 
# - SD: the separation between populations.

if (!exists("NPOPS")) NPOPS <- 5
if (!exists("SD")) SD <- 1

FUN <- function(ngenes, ncells) {
    basis <- matrix(rnorm(NPOPS * ngenes, sd=SD), ncol=NPOPS)
    p1 <- sample(NPOPS, ncells, replace=TRUE)
    p2 <- sample(NPOPS, ncells, replace=TRUE)

    src <- basis[,p1,drop=FALSE] 
    dvec <- basis[,p2,drop=FALSE] - src
    t(t(dvec) * runif(ncells)) + src
}

source("functions.R")
dir.create("results", showWarnings=FALSE)
runSimulation(sprintf("results/trajectory-%s-%s.rds", NPOPS, SD), FUN)
