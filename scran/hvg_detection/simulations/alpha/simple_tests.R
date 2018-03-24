################################################
# Simple simulations to demonstrate the accuracy of the underlying approximations
# of the distribution of the sample variances or CV2 values. 

library(matrixStats)
ngenes <- 10000
niters <- 10

thresholds <- c(0.005, 0.01, 0.05)
assessAlpha <- function(p) {
    findInterval(thresholds, sort(p))/length(p)
}

set.seed(20000)

################################################

for (mode in c("ScaledPoisson", 
               "NegativeBinomial",
               "ZeroInflatedNB")) {
    host.file <- paste0(mode, ".tsv")
    overwrite <- TRUE 

    for (ncells in c(100, 200, 500)) {       
        if (mode=="ScaledPoisson") {
            settings <- expand.grid(lambda=c(1, 10, 100, 100), scale=c(1, 2, 5))
            FUN <- function(par) { 
                par["scale"]*matrix(rpois(ngenes * ncells, lambda=par["lambda"]), ncol=ncells)
            }
        } else if (mode=="NegativeBinomial") {
            settings <- expand.grid(mu=c(1, 10, 100, 100), disp=c(0.1, 0.5, 1, 2, 5))
            FUN <- function(par) { 
                matrix(rnbinom(ngenes * ncells, mu=par["mu"], size=1/par["disp"]), ncol=ncells)
            }
        } else {
            settings <- expand.grid(mu=c(1, 10, 100, 100), zero=c(0.1, 0.2, 0.5, 0.8))
            FUN <- function(par) { 
                out <- matrix(rnbinom(ngenes * ncells, mu=par["mu"], size=10), ncol=ncells)
                out[rbinom(ngenes*ncells, 1, par["zero"])==1] <- 0L
                return(out)
            }
        }

        for (sc in seq_len(nrow(settings))) {
            current <- unlist(settings[sc,])
            collected <- vector("list", niters)

            for (it in seq_len(niters)) {
                X <- FUN(current)
                lX <- log2(X+1)
                is.spike <- 1:50 # Assuming the first 50 genes are spike-in transcripts.

                # Computing p-value by Chi-squared approximation of the distribution of CV2 estimates.
                cv <- rowVars(X)/rowMeans(X)^2
                mspike.cv <- mean(cv[is.spike])
                p.cv <- pchisq(cv/mspike.cv * (ncells-1), df=ncells-1, lower.tail=FALSE)
                alpha.cv <- assessAlpha(p.cv)
   
                # Computing p-value by Chi-squared approximation of the distribution of variance estimates of log-counts. 
                lv <- rowVars(lX)
                mspike.lv <- mean(lv[is.spike])
                p.lv <- pchisq(lv/mspike.lv * (ncells-1), df=ncells-1, lower.tail=FALSE)
                alpha.lv <- assessAlpha(p.lv)

                # Computing p-value by normality approximation of the variance-of-logs (comes out as a t-test).
                lv.spike <- lv[is.spike]
#                p.n <- sapply(lv, FUN=function(x) t.test(x, lv.spike, var.equal=TRUE)$p.value)
                n.spike <- length(lv.spike)
                p.n <- pt((lv - mspike.lv)/(sd(lv.spike)*sqrt(1+1/n.spike)), df=n.spike-1, lower.tail=FALSE)
                alpha.n <- assessAlpha(p.n)
                
                collected[[it]] <- c(current,
                                     setNames(alpha.cv, paste0("CV2_", thresholds)),
                                     setNames(alpha.lv, paste0("LogVar_", thresholds)),
                                     setNames(alpha.n, paste0("NormVar_", thresholds)))
            }

            output <- data.frame(Setting=sc, N=ncells, do.call(rbind, collected)) 
            write.table(file=host.file, output, sep="\t", quote=FALSE, row.names=FALSE, 
                        col.names=overwrite, append=!overwrite)
            overwrite <- FALSE
        }
    }
}

################################################
# Making pretty pictures.

nmethods <- 3
plotScenarios <- function(incoming, ncells) {
    blah <- read.table(incoming, header=TRUE)
    blah <- blah[blah$N==ncells,]

    by.scenario <- split(blah[,tail(seq_len(ncol(blah)), nmethods*length(thresholds))], blah$Setting)
    all.means <- lapply(by.scenario, colMeans)
    all.means <- do.call(rbind, all.means)
    all.se <- lapply(by.scenario, FUN=function(x) { sqrt(colVars(as.matrix(x))/nrow(x)) })
    all.se <- do.call(rbind, all.se)

    for (it in seq_along(thresholds)) {
        to.use <- (seq_len(nmethods) - 1)*length(thresholds) + it
        heights <- t(all.means[,to.use])
        x <- barplot(heights, ylab="Type I error rate", beside=TRUE, xlab="Simulation scenario")
        abline(h=thresholds[it], col="red", lwd=2, lty=2)

        se <- t(all.se[,to.use])
        upper <- heights + se
        segments(x, heights, x, upper)
        segments(x-0.2, upper, x+0.2, upper)

        legend("topleft", fill=grey.colors(3), rownames(heights))
    }
}

for (ncells in c(100, 200, 500)) {      
    pdf(sprintf("ScaledPoisson%i.pdf", ncells), width=12)
    plotScenarios("ScaledPoisson.tsv", ncells=ncells)
    dev.off()
}

for (ncells in c(100, 200, 500)) {      
    pdf(sprintf("NegativeBinomial%i.pdf", ncells), width=12)
    plotScenarios("NegativeBinomial.tsv", ncells=ncells)
    dev.off()
}

for (ncells in c(100, 200, 500)) {     
    pdf(sprintf("ZeroInflatedNB%i.pdf", ncells), width=12)
    plotScenarios("ZeroInflatedNB.tsv", ncells=ncells)
    dev.off()
}
