args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("<plasmidpath> <directory_for_data> <med_cov_table> \n", call.=FALSE)
} else {
  cat("Plasmid path", args[1], "\n")
}

library(ggplot2)
library(gridExtra)

plasmidpath = read.table(args[1], sep="\t",header=T)
plasmids = (unique(plasmidpath$Plasmid))
plas_count = 0
path = args[2]
data_for_names = read.table(args[3],header=T,sep='\t')
subjects = colnames(data_for_names)[-length(data_for_names)]
plot_vector = NULL
vector_count = 0
for (plas in plasmids){
  plas_count = substr(plas,2,nchar(plas)+1)
  cat(plas_count,'\n')
  uniq_plas_path = plasmidpath[which(plasmidpath$Plasmid==(plas)),]
  cat('Making graph for',plas,'\n')
  count = 1 #remove
  for (subject in subjects){
    cat(subject,'\n')
    path_to_table = paste(args[2],'/',subject,'/cov_table',sep='')
    cat(path_to_table,'\n')
    covtable = read.table(path_to_table, sep='\t',header=T)
    plasmid_cov = NULL
    plasmid_cov = covtable[which(covtable$Plasmid==(plas_count)),,]
    df = NULL
    if(count>1){
      df = covtable[which(covtable$Plasmid==(plas_count)),,]
    }
    count
    if(count == 1){
      print(head(plasmid_cov))
      plot = ggplot(data = plasmid_cov, aes(x = Base_Pair, y = Cov,size=0.1)) +
        geom_line(aes(y = Cov), size=0.05, color="darksalmon") + 
        geom_line(aes(y = Out_cyc_pos), size=0.05, color="slateblue") + 
        geom_line(aes(y = Out_cyc_neg), size=0.05, color="forestgreen") +
        theme_void()
    } else {
      print(head(df))
      plot = plot + 
        geom_line(data = df,aes(y = Cov), size=0.05, color="darksalmon") + 
        geom_line(data = df,aes(y = Out_cyc_pos), size=0.05, color="slateblue") + 
        geom_line(data = df,aes(y = Out_cyc_neg), size=0.05, color="forestgreen")
    }
    count = count + 1
  }
  vector_count = vector_count + 1
  plot = plot + ggtitle(plas_count) + theme_grey(base_size = 1.3)
  plot_vector[[vector_count]] = plot 
}
print('hereee')
final_plot = grid.arrange(grobs=plot_vector, ncol=14,nrow=14)
ggsave(filename=paste('all_good_cycs.pdf',sep=''), plot=final_plot, units="in", path=path, width=8, height=4, dpi=300)
# ggsave(filename=paste('all.pdf',sep=''), plot=last_plot(), units="in", path=path, width=8, height=4, dpi=300)
