# This script generates gaussian cluster simulations with variable gene-wise variances.
# The aim being to test whether or not Gavish and Donoho's optimal choice works when sigma is not constant.
# It expects:
# - NPOPS: the number of populations.
# - SD: the standard deviation of the population means.

if (!exists("NPOPS")) NPOPS <- 5
if (!exists("SD")) SD <- 1

FUN <- function(ngenes, ncells) {
    in.pop <- sample(NPOPS, ncells, replace=TRUE)
    pop.mean <- matrix(rnorm(NPOPS * ngenes, sd=SD), ncol=NPOPS)
    pop.mean[,in.pop,drop=FALSE]
}

source("functions.R")
dir.create("results", showWarnings=FALSE)
runSimulation(sprintf("results/genevar-%s-%s.rds", NPOPS, SD), 
    truth.FUN=FUN, 
    iters=10,
    observed.FUN=function(truth) {
        sd <- sqrt(rgamma(nrow(truth), 2, 2)) 
        truth + rnorm(length(truth), mean=truth, sd=sd)  
    }
)
