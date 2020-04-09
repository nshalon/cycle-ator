import argparse
from collections import defaultdict

#############################################################################
# Purpose: given a list of k-mer mappings, and mapping files from original megahit
# contig names to renamed contig names, output the contig coverages in matrix and list form
# to match igraph parsing (first go over rows, then columns in matrix)
# Intput: contig name mappings, kmer read mappings, adjacency list
# Output: text file with edge name and coverage, matrix with each node and coverage
#############################################################################

parser = argparse.ArgumentParser(description='rename reads to have FW and RV')
parser.add_argument('ifn_contig_names', metavar='<name map table for all original contig names>', type=str, help='input fastq files')
parser.add_argument('ifn_kmer_names', metavar='<file to connect kmer names to edges>', type=str, help='input fastq files')
parser.add_argument('ifn_kmer_map', metavar='<table of mapped k-reads to kmers>', type=str, help='input fastq files')
parser.add_argument('ifn_adj_matrix', metavar='<adjacency matrix for edge>', type=str, help='input fastq files')
parser.add_argument('ofn_edge_cov_list', metavar='<list of edge coverages>', type=str
, help='input k integer used to assemble genome')
parser.add_argument('ofn_edge_cov_matrix', metavar='<matrix of edges with values as weight>', type=str, help='fasta with the k-1mers renamed to map against ')

args = parser.parse_args()

#put all the kmers attached to their variable edges in htable
#initiate coverage of kmers to 0
edge_to_k = {}
k_coverage = {}
with open(args.ifn_kmer_names) as contig_edges:
    for line in contig_edges:
        line = line.rstrip().split(":")
        k_label = line[0]
        parent,child = line[1],line[2]
        edge_to_k[(parent,child)] = k_label
        k_coverage[k_label] = 0

#parse coverage from original megahit name
contig_to_coverage = {}
with open(args.ifn_contig_names) as contig_names:
    for line in contig_names:
        line = line.split()
        megahit_label = line[0]
        rename_contig = line[1]
        coverage = float(megahit_label.split('_')[5])
        contig_to_coverage[rename_contig] = coverage

#get coverage of variable edges
with open(args.ifn_kmer_map) as kmer_map:
    next(kmer_map)
    for line in kmer_map:
        kmer = line.split()[1]
        original_coverage = k_coverage[kmer]
        k_coverage[kmer] = original_coverage + 1

#get all nodes in alphabetical order from header of adj matrix
nodes = []
omatrix = open(args.ofn_edge_cov_matrix,"w+")
with open(args.ifn_adj_matrix) as adj_matrix:
    omatrix.write("\t")
    for line in adj_matrix:
        line = line.split()
        for node in line:
            nodes.append(node)
            omatrix.write(node + "\t")
        omatrix.write("\n")
        break
        
print(contig_to_coverage.items())
ocov = open(args.ofn_edge_cov_list,"w+")
for parent in nodes:
    omatrix.write(parent + "\t")
    for child in nodes: 
        parent_dir = parent[-1]
        child_dir = child[-1]
        if parent_dir != child_dir:
            parent_contig = parent[:3]
            child_contig = child[:3]
            if parent_contig == child_contig:
                if parent_dir == "-":
                    contig_coverage = str(contig_to_coverage[parent_contig[:2]])
                    omatrix.write(contig_coverage + "\t")
                    ocov.write(parent + " " + child + " " + contig_coverage + "\n")
                else:
                    omatrix.write("0\t")
            else:
                edge_k = edge_to_k.get((parent,child))
                if edge_k != None:
                    edge_coverage = str(k_coverage[edge_k])
                    omatrix.write(edge_coverage + "\t")
                    ocov.write(parent + " " + child + " " + edge_coverage + "\n")
                else:
                    omatrix.write("0\t")
        else:
            omatrix.write("0\t")
    omatrix.write("\n")
                





