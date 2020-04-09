args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop("<cov_table> <tail_width> <plasmid_paths> <out_dir> <stats_table> <k> \n", call.=FALSE)
} else {
  cat("Reading coverage table", args[1], "\n")
  cat("Tail width", args[2], "\n")
  cat("Reading plasmid path", args[3], "\n")
  cat("Saving to", args[4], "\n")
  cat("Reading stats table", args[5], "\n")
}

library(ggplot2)

k = as.integer(args[6])
path = args[4]
table = read.table(args[1],sep="\t",header=T)
plasmids = as.integer(unique(table$Plasmid))
plasmidpath = read.table(args[3], sep="\t",header=T)
stats_table = read.table(args[5], sep="\t",header=T)
print(plasmids)

for (i in plasmids){
  plas = paste("P",i,sep='')
  print(plas)
  uniq_plas_path = plasmidpath[which(plasmidpath$Plasmid==(plas)),]
  total_length = -k*length(uniq_plas_path[,4])
  total_length_cycle = sum(as.integer(uniq_plas_path[,4]),total_length)
  labels = uniq_plas_path[,3]
  vert_line_coords = append(uniq_plas_path[,5],(total_length_cycle))
  label_coords = vert_line_coords[-length(vert_line_coords)]+0.5*diff(vert_line_coords)
  plasmid_cov = table[which(table$Plasmid==(i)),,]
  stats = stats_table[which(stats_table$PCE==plas),]
  stat_labels = c("Z=-2","Z=-3","Z=-4","Z=-5","Z=-6","Z=-7","Half Med", "2 reads")
  plot = ggplot(data = plasmid_cov, aes(x = Base_Pair, y = Cov)) + 
    geom_line(aes(y = Cov), size=0.2, color="darksalmon") + 
    geom_line(aes(y = Out_cyc_pos), size=0.2, color="slateblue") + 
    geom_line(aes(y = Out_cyc_neg), size=0.2, color="forestgreen") + 
    ggtitle(paste(plas,' Coverage (tail width = ', args[2],")", sep='')) + 
    ylim(-1,NA) + 
    geom_vline(xintercept = vert_line_coords, colour="black", linetype = "dotted", alpha = 0.2) +
    geom_hline(yintercept = stats$Two_under, colour="black", linetype = "dotted", alpha = 0.2) +
    geom_hline(yintercept = stats$Three_under, colour="black", linetype = "dotted", alpha = 0.2) +
    geom_hline(yintercept = stats$Four_under, colour="black", linetype = "dotted", alpha = 0.2) +
    geom_hline(yintercept = stats$Five_under, colour="black", linetype = "dotted", alpha = 0.2) +
    geom_hline(yintercept = stats$Six_under, colour="black", linetype = "dotted", alpha = 0.2) +
    geom_hline(yintercept = stats$Seven_under, colour="black", linetype = "dotted", alpha = 0.2) +
    geom_hline(yintercept = 0.5*as.integer(stats$Med_cov), colour="black", linetype = "dotted", alpha = 0.2) + 
    geom_hline(yintercept = 2, colour="black", linetype = "dotted", alpha = 0.2)
  ymax = (layer_scales(plot)$y$range$range[2])
  xmax = (layer_scales(plot)$x$range$range[2]) + 20
  plot = plot + annotate("text", x=label_coords, y=c(-0.5), label=labels, size=1.5)
  plot = plot + annotate("text", x=xmax, y=c(stats$Two_under,stats$Three_under,stats$Four_under,stats$Five_under,stats$Six_under,stats$Seven_under,0.5*as.integer(stats$Med_cov),2), label=stat_labels, size=1.5)
  ggsave(filename=paste(plas,' cov.pdf',sep=''), plot=last_plot(), units="in", path=path, width=8, height=4, dpi=300)
}
