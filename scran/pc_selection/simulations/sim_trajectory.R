# This script generates trajectories between cell populations.
# It expects:
# - NPOPS: the number of populations. 
# - SD: the standard deviation of the population means.
# - HYPERVAR: a string specifying the variability of the technical noise across genes.

if (!exists("NPOPS")) NPOPS <- 5
if (!exists("SD")) SD <- 1
if (!exists("HYPERVAR")) HYPERVAR <- "none"

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
runSimulation(sprintf("results/trajectory-%i-%s-%s", NPOPS, SD, HYPERVAR), 
    truth.FUN=FUN, 
    observed.FUN=addNoise(HYPERVAR))
