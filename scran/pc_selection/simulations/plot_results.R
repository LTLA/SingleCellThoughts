# Runs through the *.rds files in 'results/' and creates pictures.
# Each plot represents the average results of one simulation scenario.

pch <- c(denoised=1, parallel=2, marchenko=5, gavish=6)
col <- c(denoised="red", parallel="blue", marchenko="goldenrod", gavish="forestgreen", optimal="grey")

dir.create("pics", showWarnings=FALSE)
for (res in list.files("results", pattern="txt$")) { 
    fname <- file.path("results", res)
    rdsname <- sub("txt$", "rds", fname)
    current <- readRDS(rdsname)
    statistics <- read.table(fname, header=TRUE)

    pdf(file.path("pics", sub(".txt$", ".pdf", res)))
    collated <- list()
    for (s in seq_len(nrow(statistics))) {
        cur.scen <- statistics[s,1:3]
        title <- paste(paste0(names(cur.scen), " = ", as.vector(cur.scen)), collapse=", ")
        cur.stats <- current[[s]]
        plot(cur.stats, pch=16, xlab="PC", ylab="Average MSE", main=title)

        cur.number <- statistics[s,-(1:3)]
        FUN <- splinefun(seq_along(cur.stats), cur.stats)
        collated[[s]] <- sapply(cur.number, FUN)

        for (x in names(pch)) { 
            points(cur.number[[x]], collated[[s]][[x]], pch=pch[x], col=col[x], cex=1.2)
        }

        if (which.max(cur.stats)==length(cur.stats)) {
            position <- "topleft"
        } else {
            position <- "topright"
        }
        legend(position, pch=pch, col=col[names(pch)], legend=names(pch))
    }
   
    collated <- do.call(rbind, collated)
    collated <- as.matrix(collated)
    scenarios <- statistics[,1:3]
    rownames(collated) <- do.call(paste, c(as.list(scenarios), sep=", "))

    par(mar=c(8.1, 4.1, 4.1, 2.1))
    barplot(t(collated), beside=TRUE, las=2, ylab="MSE", main=sub(".txt", "", res), col=col)
    legend("topright", fill=col, legend=names(col))
    dev.off()
}
