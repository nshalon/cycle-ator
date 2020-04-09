args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("<summary_table> <out stats directory> <out_graph_dir> <full graph> <edge_weights> \n", call.=FALSE)
} else {
  cat("summary table", args[1], "\n")
  cat("out stats directory", args[2], "\n")
  cat("out graph directory", args[3], "\n")
}

library(igraph)
library(reshape2)
library(dplyr)
#library(ggplot2)

# odir = paste(path,"/cycle_test/out_dir/",sep="")
odir = args[2]

summary = read.delim(args[1],sep='\t',header=T)
cycle_summary = as.data.frame(summary[match(unique(summary$cycle), summary$cycle),])

full_edge_summary = as.data.frame(read.delim(args[5],sep="\t",header=T))
edge_weights = as.character(round(as.numeric(full_edge_summary[,5]),2))

all_matrix = read.table(args[4])
all_names = as.vector(rownames(all_matrix))
all_names = substring(all_names,2)
colnames(all_matrix) = rownames(all_matrix) = all_names
all_matrix = as.matrix(all_matrix)

network = graph_from_adjacency_matrix(all_matrix)

perCycle <- function(row) {
  
  rowlist = as.data.frame(row)
  cycle = rowlist["cycle",]
  adj_matrix = read.delim(paste(odir,"/",cycle,"_adjmatrix.txt",sep=""),sep="\t")
  
  mat_df = as.data.frame(adj_matrix)
  names = mat_df$X
  mat_df$X.1 = NULL
  mat_df = mat_df[,2:ncol(mat_df)]
  names = substring(names,2)
  rownames(mat_df) = colnames(mat_df) = names
  cycle_edges = names
  
  mat = as.matrix(mat_df)
  
  total_out = as.character(rowlist["total_out_cov",])
  total_in = as.character(rowlist["total_in_cov",])
  avg_cov = as.character(rowlist["avg_cov",])
  bottleneck_cov = as.character(rowlist["bottleneck_cov",])
  cyc_length = as.character(rowlist["cycle_length",])
  summary = paste("Total Out Cov:",total_out,"Total In Cov:",total_in,
              "Average Cov:",avg_cov,"Bottleneck Cov:",bottleneck_cov,
              "Cycle Length:",cyc_length)
  # edge_labs = read.delim(paste(odir,"/",cycle,"_edgeweights.txt",sep=""),sep="\t",header=F)
  # edge_labs = as.data.frame(edge_labs)
  # colnames(edge_labs) = c("edge","stat")
  # edge_labs = edge_labs[order(edge_labs$edge),]
  # edge_weights_cyc = as.character(edge_labs$stat)
  colors = c("orange","red")
  node_color_vec = colors[1 + all_names %in% names]
  file = paste(args[3],"/",cycle,"_network_graph.pdf",sep="")
  pdf(file,width=7,height=9)
  
  plot(network,
        vertex.size = 10,
        vertex.label.cex = 1,
        vertex.color = node_color_vec,
        edge.arrow.size = 0.5,
        edge.label = edge_weights,
        edge.label.cex = 1
       )
                      
  title(main = cycle, sub = summary,
        cex.sub = 0.7)
  
  dev.off()

}

apply(cycle_summary, 1, perCycle)



