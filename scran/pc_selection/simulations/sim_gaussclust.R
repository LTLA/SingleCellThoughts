# This script generates gaussian cluster simulations.
# It expects:
# - NPOPS: the number of populations.
# - SD: the standard deviation of the population means.
# - HYPERVAR: a string specifying the variability of the technical noise across genes.

if (!exists("NPOPS")) NPOPS <- 5
if (!exists("SD")) SD <- 1
if (!exists("HYPERVAR")) HYPERVAR <- "none"

FUN <- function(ngenes, ncells) {
    in.pop <- sample(NPOPS, ncells, replace=TRUE)
    pop.mean <- matrix(rnorm(NPOPS * ngenes, sd=SD), ncol=NPOPS)
    pop.mean[,in.pop,drop=FALSE]
}

source("functions.R")
dir.create("results", showWarnings=FALSE)
runSimulation(sprintf("results/gaussclust-%i-%s-%s", NPOPS, SD, HYPERVAR), 
    truth.FUN=FUN, 
    observed.FUN=addNoise(HYPERVAR))
