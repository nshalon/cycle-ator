args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("<summary_cycles> <out_dir> \n", call.=FALSE)
} else {
  cat("summary of cycles file:",args[1])
  cat("out dir", args[2], "\n")
}
path = args[2]
library(ggplot2)


cycles = read.table(args[1],sep="\t",header=T)
cycles = cycles[!duplicated(cycles$cycle),]
cycles$maxFlow = apply(cycles[,c(2,3)], 1, max)

plot = ggplot(cycles,aes(log10(maxFlow),log10(bottleneck_cov))) + 
  geom_jitter(height=1,width=1) + 
  geom_abline(slope=1,linetype="dotted") +
  ylim(0,NA) +
  xlim(0,NA) +
  ylab("Bottleneck Coverage") +
  xlab("Max of In/Out Coverage")

ggsave("cycle_coverages.pdf",plot=plot,path=path)







  

