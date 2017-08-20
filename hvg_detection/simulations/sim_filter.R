library(scran)
source("functions.R")

set.seed(200)
ngenes <- 1000
ncells <- 1000
nspikes <- 1000

pdf("results_filter.pdf")
boundaries <- c(-20, 5)
for (scale in c(0, 0.001, 0.01)) { 
    for (asymptote in c(0.1, 0.2, 0.5)) { 
        sim.spikes <- sampleCounts(ngenes=nspikes, nsamples=ncells, bounds=boundaries, scale=scale, asymptote=asymptote)
        is.spike <- seq_len(nspikes) + ngenes
        lcounts <- log2(sim.spikes+1)
        fit <- trendVar(lcounts, span=0.2)
        plot(fit$mean, fit$var/fit$trend(fit$mean), log="x", main=sprintf("Scale of %.3f, asymptote of %.3f", scale, asymptote))
        abline(v=0.1, col="red", lty=2)
    }
}
dev.off()
