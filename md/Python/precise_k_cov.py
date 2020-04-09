import argparse
from collections import defaultdict
import sys

#############################################################################
# Purpose: given a list of k-mer mappings, and mapping files from original megahit
# contig names to renamed contig names, output the contig coverages in matrix and list form
# to match igraph parsing (first go over rows, then columns in matrix)
# Intput: contig name mappings, kmer read mappings, adjacency list
# Output: text file with edge name and coverage, matrix with each node and coverage
#############################################################################

parser = argparse.ArgumentParser(description='rename reads to have FW and RV')
parser.add_argument('ifn_R1', metavar='<first fastq of paired reads>', type=str, help='input fastq files')
parser.add_argument('ifn_R2', metavar='<second fastq paired reads>', type=str, help='input fastq files')
parser.add_argument('ifn_contig_names', metavar='<name map table for all original contig names>', type=str, help='input fastq files')
parser.add_argument('ifn_kmer_names', metavar='<file to connect kmer names to edges>', type=str, help='input fastq files')
parser.add_argument('ifn_kmer_reads', metavar='<table of mapped k-reads to kmers>', type=str, help='input fastq files')
parser.add_argument('ifn_adj_list', metavar='<adjacency list for edge>', type=str, help='input fastq files')
parser.add_argument('ifn_k_fasta', metavar='<fasta of all k edges>', type=str, help='fasta that has names and sequences of all k edges')
parser.add_argument('ofn_edge_cov_list', metavar='<list of edge coverages>', type=str, help='input k integer used to assemble genome')
parser.add_argument('ofn_edge_cov_matrix', metavar='<matrix of edges with values as weight>', type=str, help='fasta with the k-1mers renamed to map against ')

args = parser.parse_args()

def reverse_complement(seq):
    comp_base = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
    comp_seq = ''
    for base in seq[::-1]:
        comp_seq += comp_base[base]
    return comp_seq

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

parent_to_child = {}
with open(args.ifn_adj_list) as adj_list:
    for line in adj_list:
        line = line.rstrip().split(":")
        if len(line[1]) <= 1: #parent node with no children
            continue
        parent,children = line[0],line[1]
        for child in children.split(","):
            parent_to_child.setdefault(parent, []).append(child)

#parse coverage from original megahit name
contig_to_coverage = {}
contig_to_length = {}
with open(args.ifn_contig_names) as contig_names:
    for line in contig_names:
        line = line.rstrip().split()
        megahit_label = line[0]
        rename_contig = line[1]
        coverage = float(megahit_label.split('_')[5])
        contig_to_coverage[rename_contig] = coverage
        contig_to_length[rename_contig] = (megahit_label.split('_')[3])

k_sequences_to_label = {}
k_size = 0

with open(args.ifn_k_fasta) as k_fasta:
    for line in k_fasta:
        line = line.rstrip()
        if line[0] == ">":
            k_label = line[1:]
        else:
            k_sequence = line
            k_sequences_to_label[k_sequence] = k_label

k_length = str(len(k_sequence))

print("kmer size:", k_length)

#get coverage of variable edges
with open(args.ifn_kmer_reads) as k_reads:
    next(k_reads)
    fastq_line_count = 0
    k_map_count = 0 
    for line in k_reads:
        if (fastq_line_count % 4000000) == 0:
            print("Kmer reads count",k_reads_count)
        if fastq_line_count % 4 == 1:
            k_reads_count += 1
            sequence = line.rstrip()
            comp_sequence = reverse_complement(sequence) 
            if k_sequences_to_label.get(sequence) != None:
                k_label = k_sequences_to_label[sequence]
                k_coverage[k_label] += 1
                k_map_count += 1
            elif k_sequences_to_label.get(comp_sequence) != None:
                k_label = k_sequences_to_label[comp_sequence]
                k_coverage[k_label] += 1
                k_map_count += 1    
        fastq_line_count += 1
        

print("Total k reads",k_reads_count)
print("Percent k reads match to kmer:", str(100*round(float(k_map_count/k_reads_count),10)) + "%")

        
ocov = open(args.ofn_edge_cov_list,"w+")
ocov.write("parent\tchild\tcoverage\tlength\tweight\n")
for parent in parent_to_child.keys():
    parent_dir = parent[-1]
    children = parent_to_child[parent]
    if parent_dir == "-":
        invariable_edge = True
    else:
        invariable_edge = False
    for child in children:
        if invariable_edge:
            parent_contig = parent[:-1]
            child_contig = child[:-1]
            contig_coverage = str(contig_to_coverage[parent_contig[:-1]])
            contig_length = contig_to_length[parent_contig[:-1]]
            contig_weight = str(round(float(contig_length) / float(contig_coverage),3))
            ocov.write(parent + "\t" + child + "\t" + contig_coverage + "\t" + contig_length + "\t" + contig_weight + "\n")
        else:
            edge_k = edge_to_k.get((parent,child))
            edge_coverage = str(k_coverage[edge_k])
            if edge_coverage == "0":
                edge_coverage = str(0.00000001)
            weight = str(round(float(k_length) / float(edge_coverage),10))
            ocov.write(parent + "\t" + child + "\t" + edge_coverage + "\t" + str(k_length) + "\t" + weight + "\n")