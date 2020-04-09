args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop("<num_subjects> <plasmidpath> <directory_for_data> <out_dir> <k> <med_cov_table> \n", call.=FALSE)
} else {
  cat("Number of subejects", args[1], "\n")
  cat("Plasmid path", args[2], "\n")
}

library(ggplot2)
library(gridExtra)


data_for_names = read.table(args[6],header=T,sep='\t')
subjects = colnames(data_for_names)[-length(data_for_names)]
k = as.integer(args[5])
num_subjects = args[1]
plasmidpath = read.table(args[2], sep="\t",header=T)
plasmids = (unique(plasmidpath$Plasmid))
print(plasmids)
plas_count = 0
path = args[4]
plot_vector = NULL
for (plas in plasmids){
  plas_count = substr(plas,2,nchar(plas)+1)
  uniq_plas_path = plasmidpath[which(plasmidpath$Plasmid==(plas)),]
  total_length = -k*length(uniq_plas_path[,4])
  total_length_cycle = sum(as.integer(uniq_plas_path[,4]),total_length)
  labels = uniq_plas_path[,3]
  vert_line_coords = append(uniq_plas_path[,5],(total_length_cycle))
  label_coords = vert_line_coords[-length(vert_line_coords)]+0.5*diff(vert_line_coords)
  plot_vector = plot_vector[F]
  cat('Making graph for',plas,'\n')
  subject_count = 0
  for (subject in subjects){
    subject_count = subject_count + 1
    path_to_table = paste(args[3],'/',subject,'/cov_table',sep='')
    covtable = read.table(path_to_table, sep="\t",header=T)
    plasmid_cov = covtable[which(covtable$Plasmid==(plas_count)),,]
    plot = ggplot(data = plasmid_cov, aes(x = Base_Pair, y = Cov, size=0.5)) + 
      geom_line(aes(y = Cov), size=0.2, color="darksalmon") + 
      geom_line(aes(y = Out_cyc_pos), size=0.2, color="slateblue") + 
      geom_line(aes(y = Out_cyc_neg), size=0.2, color="forestgreen") +
      ggtitle(paste(plas,' Coverage, ', subject ,sep='')) + 
      theme_grey(base_size = 1.3) +
      ylim(-1,NA) + 
      geom_vline(xintercept = vert_line_coords, colour="black", linetype = "dotted", alpha = 0.2)
    plot = plot + annotate("text", x=label_coords, y=c(-0.5), label=labels, size=1)
    plot_vector[[subject_count]] = plot
    cat('length vector',length(plot_vector))
  }
  final_plot = grid.arrange(grobs=plot_vector, nrow=5, ncol=4)
  # final_plot = grid.arrange(p1, p2, p3, p4, p5, ncol=1, heights=unit(c(1,1,1,1,1), c("in","in","in","in")),newpage=T)
  # print(length(plot_vector))
  ggsave(filename=paste(plas,'_cov_all_subjects.pdf',sep=''), plot=final_plot, units="in", path=path, width=8, height=4, dpi=300)
} 
