args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("<coverage table> <out plots directory> \n", call.=FALSE)
} else {
  cat("coverage table", args[1], "\n")
  cat("out plots directory", args[2], "\n")
}

contig_cov = read.table(args[1],header=T)

library(ggplot2)

plot = ggplot(contig_cov,aes(contig,coverage)) +
  geom_bar(stat="identity")


ggsave(filename="contig_coverage.pdf", plot=plot, units="in", path=args[2], width=8, height=4, dpi=300)
