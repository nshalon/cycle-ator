args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("<summary> <subjects> <directory_for_data> <out_dir> \n", call.=FALSE)
} else {
  cat('summary',args[1])
  cat('subjects',args[2])
  cat('in_dir',args[3])
  cat('out_dir',args[4])
}

library(ggplot2)
library(gridExtra)

plasmids = read.table(args[1],header=T,sep='\t')
plasmids = (plasmids[,1])
plasmids = as.character(plasmids)
plasmids = as.integer(substr(plasmids,2,nchar(plasmids)+1))
subjects = read.table(args[2],header=F,sep='\t')
subjects = as.vector(subjects[,1])
plot_vector = NULL
colors = c('red4','yellow3','green4','blue4')
for(subject in subjects){
  path = paste(args[4])
  snp_table_path = paste(args[3],'/',subject,sep='')
  snp_table = read.table(snp_table_path,header=T,sep='\t')
  plot_vector = plot_vector[F]
  plas_count = 0
  cat(paste('\nmaking graph for subject',subject))
  for(plas in plasmids){
    plas_count = plas_count + 1
    plasmid_cov = snp_table[which(snp_table$cycle==(plas)),,]
    plot = ggplot(data=plasmid_cov,aes(x=t1,y=t2)) + 
      geom_point(aes(color=base,size=),size=0.00001,alpha=0.5) + 
      ggtitle(paste(subject,plas)) + 
      theme_grey(base_size = 1.3) +
      scale_colour_manual(values = colors) + 
      theme(legend.position = "none") +
      xlim(0,1) + 
      ylim(0,1)
    plot_vector[[plas_count]] = plot
  }
  final_plot = grid.arrange(grobs=plot_vector, nrow=14, ncol=14)
  cat('past the final plot')
  ggsave(filename=paste(subject,'longitudinal snps.pdf'), plot=final_plot, units="in", path=path, width=8, height=4, dpi=300)
}

