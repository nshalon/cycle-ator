gen.read.coords=function(genome.table.ifn, insert=200, insert.sd=100, read.length=150, ofn)
{
    options(stringsAsFactors=F)
    table = read.delim(genome.table.ifn)
    contigs = table$contig

    result = NULL

    total.contig.size = sum(table$length)
    N.total = ceiling(sum(table$xcov*table$length) / (2*read.length))
    cat(sprintf("covering %.1fkb of contigs with x-coverage in the range of %.1f to %.1f, expected number of reads: %d\n",
                total.contig.size/10^3, min(table$xcov), max(table$xcov), N.total))

    M = length(contigs)
    if (M == 0)
        stop("no contigs found")

    cat(sprintf("going over %d contigs...\n", length(contigs)))
    for (i in 1:M) {
        contig = table$contig[i]
        contig.size = table$length[i]
        xcov = table$xcov[i]
        circ = table$circ[i]

        N = ceiling((xcov*contig.size) / (2*read.length))
        start = floor(runif(N, min=0, max=contig.size))
        end = pmax(start + 1, start + 2*read.length + insert + round(rnorm(N, sd=insert.sd)))
        strand = ifelse(runif(N) > 0.5, 1, -1)

        start1 = ifelse(strand ==  1, start, end-read.length)
        start2 = ifelse(strand == -1, start, end-read.length)

        df = data.frame(
            contig1=contig, start1=start1, end1=start1+read.length, strand1=strand,
            contig2=contig, start2=start2, end2=start2+read.length, strand2=-strand)

        if (!circ) {
            df = df[pmax(df$start1,df$start2,df$end1,df$end2) < contig.size,]
        }
        result = rbind(result, df)
    }
    cat(sprintf("saving %d reads to file: %s\n", dim(result)[1], ofn))
    write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

utest.generate.read.coords.table=function()
{
    gen.read.coords(genome.table.ifn="x.tab", base.length=100, insert=300, insert.sd=200, read.length=150, ofn="coord.tab")
}
