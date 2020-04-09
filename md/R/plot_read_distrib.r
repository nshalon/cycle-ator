fig.dir=function(dir, verbose=T)
{
    if (verbose) cat(sprintf("figure dir: %s\n", dir))
    if (!file.exists(dir)) {
        command = paste("mkdir -p", dir)
        if (system(command) != 0)
            stop(sprintf("failed command: %s\n"))
    }
}

fig.start=function(ofn, type="pdf", fdir=NA, verbose=T, width=4, height=4, ...)
{
    if (!is.na(fdir))
        fig.dir(fdir)

    if (verbose) cat(sprintf("creating figure: %s\n", ofn))
    switch(type,
           png = png(ofn, width=width, height=height, ...),
           pdf = pdf(ofn, width=width, height=height, ...)
           )
}

fig.end=function()
{
    dev.off()
}

plot.read.distrib=function(ifn, fdir)
{
    df =  read.delim(ifn)
    if (!all(df$contig1 == df$contig2))
        stop("inter-contig reads")

    df$contig = df$contig1
    df$start = pmin(df$start1, df$start2)
    df$end = pmax(df$end1, df$end2)
    df$length = df$end - df$start
    contigs = unique(df$contig)

    # plot density
    ofn = paste(fdir, "/read_size.pdf", sep="")
    fig.start(ofn, fdir=fdir, width=4, height=4)
    plot(density(df$length), main="molecule length distribution")
    fig.end()

    for (contig in contigs) {
        ofn = paste(fdir, "/", contig, ".pdf", sep="")
        fig.start(ofn, width=4, height=4)
        start = df$start[df$contig == contig]
        end = df$end[df$contig == contig]
        hist(start, 100, border=NA, col="darkblue")
        fig.end()
    }
}
