# Runs through the *.rds files in 'results/' and creates pictures.
# Each plot represents the average results of one simulation scenario.
# Color encodes the MSE (yellow = lower MSE), while the size encodes the number of PCs (smaller = fewer).

dir.create("pics", showWarnings=FALSE)
library(viridis)

for (res in list.files("results")) { 
    fname <- file.path("results", res)
    current <- readRDS(fname)
    pdf(file.path("pics", sub(".rds$", ".pdf", res)))

    for (s in seq_along(current$scenarios)) {
        cur.scen <- current$scenarios[[s]]
        title <- paste(paste0(names(cur.scen), " = ", as.vector(cur.scen)), collapse=", ")

        cur.stats <- current$statistics[[s]]
        by.segment <- cut(-cur.stats$MSE, 100)
        plot(cur.stats$Bias, cur.stats$Variance, pch=16,
            xlab="Approximate bias", ylab="Estimated variance", main=title,
            col=viridis(100)[by.segment],  
            cex=sqrt(seq_along(by.segment))/5 + 1)

        cur.number <- current$retained[[s]]
        BIASFUN <- splinefun(seq_along(cur.stats$Bias), cur.stats$Bias)
        VARFUN <- splinefun(seq_along(cur.stats$Variance), cur.stats$Variance)
        points(BIASFUN(cur.number$denoised), VARFUN(cur.number$denoised), col="red", cex=1.2)
        points(BIASFUN(cur.number$parallel), VARFUN(cur.number$parallel), col="red", pch=2, cex=1.2)
        points(BIASFUN(cur.number$upper), VARFUN(cur.number$upper), col="red", pch=5, cex=1.2)
        points(BIASFUN(cur.number$optimal), VARFUN(cur.number$optimal), col="red", pch=6, cex=1.2)

        if (which.max(cur.stats$Bias)==which.max(cur.stats$Variance)) { 
            position <- "topleft"
        } else {
            position <- "topright"
        }
        legend(position, pch=c(1,2,5,6), legend=names(cur.number), col="red")
    }
    dev.off()
}
