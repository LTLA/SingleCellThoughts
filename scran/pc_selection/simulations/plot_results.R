# Runs through the *.rds files in 'results/' and creates pictures.
# Each plot represents the average results of one simulation scenario.

dir.create("pics", showWarnings=FALSE)
for (res in list.files("results", pattern="txt$")) { 
    fname <- file.path("results", res)
    rdsname <- sub("txt$", "rds", fname)
    current <- readRDS(rdsname)
    statistics <- read.table(fname, header=TRUE)

    pdf(file.path("pics", sub(".txt$", ".pdf", res)))
    for (s in seq_len(nrow(statistics))) {
        cur.scen <- statistics[s,1:3]
        title <- paste(paste0(names(cur.scen), " = ", as.vector(cur.scen)), collapse=", ")
        cur.stats <- current[[s]]
        plot(cur.stats, pch=16, xlab="PC", ylab="Average MSE", main=title)

        cur.number <- statistics[s,-(1:3)]
        FUN <- splinefun(seq_along(cur.stats), cur.stats)
        points(cur.number$denoised, FUN(cur.number$denoised), col="red", cex=1.2)
        points(cur.number$parallel, FUN(cur.number$parallel), col="red", pch=2, cex=1.2)
        points(cur.number$marchenko, FUN(cur.number$marchenko), col="red", pch=5, cex=1.2)
        points(cur.number$gavish, FUN(cur.number$gavish), col="red", pch=6, cex=1.2)

        if (which.max(cur.stats)==length(cur.stats)) {
            position <- "topleft"
        } else {
            position <- "topright"
        }
        legend(position, pch=c(1,2,5,6), col="red", 
            legend=c("denoised", "parallel", "marchenko", "gavish"))
    }
    dev.off()
}
