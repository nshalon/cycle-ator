# #!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("<class_table> <plasmid_paths> <bucket> <out_dir> \n", call.=FALSE)
} else {
  cat("Reading table", args[1], "\n")
  cat("Reading plasmid path", args[2], "\n")
  cat("Bucket", args[3], "\n")
  cat("Saving to", args[4], "\n")
}

library(ggplot2)

bucket = as.integer(args[3])
table = read.table(args[1],sep="\t",header=T)
plasmidpath = read.table(args[2], sep="\t",header=T)
path = args[4]

print(paste("Path:", path))
print("Making Graphs...")
plasmids = as.integer(unique(table$Plasmid))
plaslengths = NULL
for (i in plasmids){
  plas = paste("P",i,sep="")
  print(plas)
  uniq_plas_path = plasmidpath[which(plasmidpath$Plasmid==(plas)),]
  labels = uniq_plas_path[,3]
  for (j in c(1:2)){
    if (j==1){
      sign = '+'
    }
    else{
      sign = '-'
    }
    data = table[which(table$Plasmid==(i) & table$Strand==sign),]
    cat(length(data$Base_Pair))
    plaslengths = append(plaslengths, tail(data$Base_Pair, n=1))
    print(tail(data$Base_Pair, n=1))
    vert_line_coords = append(uniq_plas_path[,5],tail(data$Base_Pair, n=1))
    label_coords = vert_line_coords[-length(vert_line_coords)]+0.5*diff(vert_line_coords)
    title = paste("Plasmid ", toString(i), " Strand ", sign, " Bin ", bucket, sep='')
    plot = ggplot(data = data, aes(x = Base_Pair)) +
      geom_line(aes(y = Good_Reads/bucket), color = "blue", size=0.1, alpha=0.35)+
      geom_line(aes(y = Distant_Reads/bucket), color = "red", size=0.1, alpha=0.8) +
      geom_line(aes(y = Not_on_plasmid/bucket), color = "green", size=0.1, alpha=0.8) +
      geom_line(aes(y = On_lin_seq/bucket), color = "black",  size=0.1, alpha=0.8) +
      xlab('Position (bp)') +
      ylab('Read Density') + 
      ggtitle(title) +
      geom_vline(xintercept = vert_line_coords, colour="black", linetype = "dotted", alpha = 0.2)
    ymax = (layer_scales(plot)$y$range$range[2])
    plot = plot + annotate("text", x=label_coords, y=-0.04*ymax, label = labels, size=1.5)
    ggsave(filename=paste(title,'.pdf',sep=''), plot=last_plot(), units="in", path=path, width=8, height=4, dpi=300)
    cat("Graph", i, sign, "made!\n")
  }
}

plasmids = paste("P",plasmids,sep='')
plaslengths = plaslengths[c(TRUE,FALSE)]
df = data.frame("Plasmids" = plasmids, "Lengths" = plaslengths)
ggplot(data=df, aes(x=Plasmids,y=Lengths)) + geom_bar(stat = "identity",width=0.65) + ggtitle("Summary of Plasmid Lengths")
ggsave(filename="Summary.pdf", plot=last_plot(), units="in", path=path, width=8, height=4, dpi=300)

