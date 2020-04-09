args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop("<plasmidpath> <out_dir> <subject_file> <base_out> <summary> <base_out_cov> \n", call.=FALSE)
} else {
  cat("Number of subejects", args[1], "\n")
  cat("Plasmid path", args[2], "\n")
}

library(ggplot2)
library(gridExtra)

plasmidpath = read.table(args[1], sep="\t",header=T)
data_for_names = read.table(args[3],header=T,sep='\t')
subjects = colnames(data_for_names)[-length(data_for_names)]
# plasmids = (unique(plasmidpath$Plasmid))
plas_count = 0
path = args[2]
plot_vector = NULL
plasmids = (read.table(args[5],header=T,sep='\t'))$Plasmid
plasmid_lengths = (read.table(args[5],header=T,sep='\t'))$Length
print(plasmids)
print(length(plasmids))
plas_counter = 0
for (plas in plasmids){
  plas_counter = plas_counter + 1
  plas_count = substr(plas,2,nchar(plas)+1)
  print(plas_count)
  plot_vector = plot_vector[F]
  cat('Making graph for',plas,'\n')
  subject_count = 0
  as.integer(plas_count)
  plasmid_length = plasmid_lengths[plas_counter]
  cat('plas length,',plasmid_length,'\n')
  for (subject in subjects){
    path_to_table = paste(args[6],'/',subject,'/cov_table',sep='')
    covtable = read.table(path_to_table, sep="\t",header=T)
    plasmid_cov = covtable[which(covtable$Plasmid==(plas_count)),,]
    subject_count = subject_count + 1
    path_to_snp_table = paste(args[4],'/',subject,'/parsed_snp.tab',sep='')
    full_snp_table = read.table(path_to_snp_table,sep='\t',header=T,row.names=NULL)
    colnames(full_snp_table) = c("Cycle","coord","A","C","G","t")
    snp_table_for_plas = full_snp_table[which(full_snp_table$Cycle==plas_count),,]
    df_for_snp = as.data.frame(t(as.data.frame(rbind(x=c(1:plasmid_length),y=c(1:plasmid_length)))))
    plot = ggplot(data = plasmid_cov, aes(x=Base_Pair,y=Cov)) +
      geom_line(size=0.2,color='grey60') +
      geom_point(data=snp_table_for_plas,aes(x=coord,y=A), size=0.1, alpha=0.5,color="red4") + 
      geom_point(data=snp_table_for_plas,aes(x=coord,y=t), size=0.1, alpha=0.5,color="blue4") + 
      geom_point(data=snp_table_for_plas,aes(x=coord,y=G), size=0.1, alpha=0.5,color="green4") + 
      geom_point(data=snp_table_for_plas,aes(x=coord,y=C), size=0.1, alpha=0.5,color="yellow3") + 
      ggtitle(paste(plas,' SNP, ', subject ,sep='')) +
      ylab('Base_count') +
      xlab('BP') +
      theme_grey(base_size = 2) +
      xlim(0,plasmid_length)
    ymax = (ggplot_build(plot)$layout$panel_params[[1]]$y.range)[2] + 1
    plot = plot+ylim(0.1,ymax)
    plot_vector[[subject_count]] = plot
  }
  final_plot = grid.arrange(grobs=plot_vector, nrow=5, ncol=4)
  ggsave(filename=paste(plas,'_snp.pdf',sep=''), plot=final_plot, units="in", path=path, width=8, height=4, dpi=300)
} 
