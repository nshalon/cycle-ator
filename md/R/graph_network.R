args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("<adjacency matrix> <edge weights> <out plots directory> \n", call.=FALSE)
} else {
  cat("adjacency matrix", args[1], "\n")
  cat("edge weights", args[2], "\n")
  cat("out plot directory", args[3], "\n")
}

library(igraph)
library(reshape2)
library(ggplot2)

path = args[3]
# matrix = read.table("/Users/nitan/Relman/new/test_files/adjacency_matrix")
# list_cov = read.delim("/Users/nitan/Relman/new/test_files/edge_summary",sep=c("\t"," "),header=T)
list_cov = as.data.frame(read.delim(args[2],sep="\t",header=T))
edge_labs = paste0("l:",list_cov$length,",c:",round(as.numeric(list_cov$coverage)),",w:",signif(as.numeric(list_cov$weight),digits=1))
# edge_weights = as.character(round(as.numeric(list_cov[,5]),2))

matrix = read.table(args[1])
names = as.vector(rownames(matrix))
names = substring(names,2)
colnames(matrix) = rownames(matrix) = names
matrix = as.matrix(matrix)

network = graph_from_adjacency_matrix(matrix)


pdf(paste(path,"/network_graph.pdf",sep=""),width=7,height=9)

network_plot = plot(network,
     vertex.size=20,
     vertex.label.cex = 0.8,
     edge.arrow.size=0.5,
     edge.label = edge_labs,
     edge.label.cex=0.35)
dev.off()

tile = melt(matrix)
matrix_plot = ggplot(tile, aes(x=Var2,y=Var1,fill=value)) + 
  geom_tile() +
  scale_y_discrete(limits = rev(levels(tile$Var2))) + 
  scale_x_discrete(position = "top") + 
  ylab(NULL) +
  xlab(NULL) + 
  theme(legend.position = "none")

ggsave(filename="adjacency_matrix.pdf", plot=matrix_plot, units="in", path=path, width=8, height=8, dpi=300)
