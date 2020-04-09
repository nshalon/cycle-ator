args = commandArgs(trailingOnly=TRUE)
if (length(args)!=8) {
  stop("<plasmidpath> <directory_for_data> <out_dir> <k> <med_cov_table> <gene_table> <gene_table_hmm> <gc_file> \n", call.=FALSE)
} else {
  cat("Number of subejects", args[1], "\n")
  cat("Plasmid path", args[2], "\n")
}

library(ggplot2)
library(gridExtra)

gc_table = read.table(args[8],header=T,sep='\t')
hmm_table = read.table(args[7],header=T,sep=',')
gene_table = read.table(args[6],header=T,sep=',')
starts = as.integer(gene_table$start_coord)
stops = as.integer(gene_table$stop)
labels = gene_table$prot_desc
label_coords = (starts+stops)/2
data_for_names = read.table(args[5],header=T,sep='\t')
subjects = colnames(data_for_names)[-length(data_for_names)]
k = as.integer(args[4])
plasmidpath = read.table(args[1], sep="\t",header=T)
plasmids = (unique(plasmidpath$Plasmid))
print(plasmids)
path = args[3]
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
  hmm_genes = NULL
  hmm_gene_labels = NULL
  hmm_genes = hmm_table[which(hmm_table$PCE_num==plas_count),]
  hmm_starts = as.integer(hmm_genes$start_coord)
  hmm_stops = as.integer(hmm_genes$stop_coord)
  hmm_gene_labels = hmm_genes$prot_desc
  hmm_gene_label_coords = (hmm_starts+hmm_stops)/2
  hmm_add_gene = F
  if (length(hmm_gene_labels)!=0){
    hmm_add_gene = T
  }
  genes = gene_table[which(gene_table$PCE_num==plas_count),]
  starts = as.integer(genes$start_coord)
  stops = as.integer(genes$stop_coord)
  gene_labels = genes$prot_desc
  gene_label_coords = (starts+stops)/2
  add_gene = F
  if (length(gene_labels)!=0){
    add_gene = T
  }
  for (subject in subjects){
    gc_in = gc_table[which(gc_table$cycle==plas_count),]
    subject_count = subject_count + 1
    path_to_table = paste(args[2],'/',subject,'/cov_table',sep='')
    covtable = read.table(path_to_table, sep="\t",header=T)
    plasmid_cov = covtable[which(covtable$Plasmid==(plas_count)),,]
    print(subject)
    y_total = max(plasmid_cov$Cov)
    
    plot = ggplot(data = plasmid_cov, aes(x = Base_Pair, y = Cov, size=0.5)) + 
      geom_line(aes(y = Cov), size=0.2, color="darksalmon") + 
      geom_line(aes(y = Out_cyc_pos), size=0.2, color="slateblue") + 
      geom_line(aes(y = Out_cyc_neg), size=0.2, color="forestgreen") +
      ggtitle(paste(plas,' Coverage, ', subject ,sep='')) + 
      theme_grey(base_size = 1.3) +
      ylim(-1.2,NA) + 
      geom_vline(xintercept = vert_line_coords, colour="black", linetype = "dotted", alpha = 0.2) 
    if (add_gene == T){
      plot = plot + geom_segment(data=genes,aes(x=start_coord, y=c(-0.5), xend=stop_coord, yend = c(-0.5)),size=0.1) +
        geom_point(data=genes,aes(x=start_coord,y=c(-0.5)),colour='red',size=0.1) + 
        annotate(geom="text", x=gene_label_coords, y=c(-0.7), label=gene_labels,size=0.3)
    }
    if (hmm_add_gene == T){
      print('hmm')
      plot = plot + geom_segment(data=hmm_genes,aes(x=start_coord, y=c(-0.9), xend=stop_coord, yend = c(-0.9)),size=0.1) +
        geom_point(data=hmm_genes,aes(x=start_coord,y=c(-0.9)),colour='blue',size=0.1) + 
        annotate(geom="text", x=hmm_gene_label_coords, y=c(-1.1), label=hmm_gene_labels,size=0.3)
    }
    if(y_total==0){
      y_total = 1
    }
    gc_in$gc = (gc_in$gc)*y_total
    print(y_total*gc_in$gc[1:10])
    plot = plot + 
      geom_line(data=gc_in,aes(x=bp,y=(gc)),colour='turquoise3',size=0.2,alpha=0.5) +
      scale_y_continuous(sec.axis = sec_axis((~ . / y_total),name = "GC content"))
    plot_vector[[subject_count]] = plot
  }
  cat('\nhello')
  final_plot = grid.arrange(grobs=plot_vector, nrow=5, ncol=4)
  ggsave(filename=paste(plas,'.pdf',sep=''), plot=final_plot, units="in", path=path, width=8, height=4, dpi=300)
} 
