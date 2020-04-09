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
parser.add_argument('ifn_adj_matrix', metavar='<adjacency list for edge>', type=str, help='input fastq files')
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

connected = {}
with open(args.ifn_adj_list) as adj_list:
    for line in adj_list:
        line = line.rstrip().split(":")
        if len(line[1]) <= 1: #parent node with no children
            continue
        parent,children = line[0],line[1]
        for child in children.split(","):
            connected[(parent,child)] = True

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
with open(args.ifn_k_fasta) as k_fasta:
    for line in k_fasta:
        line = line.rstrip()
        if line[0] == ">":
            k_label = line[1:]
        else:
            k_sequence = line
            k_sequences_to_label[k_sequence] = k_label

k_length = len(k_sequence)

print("kmer size:", k_length)

print("Now going through R1...")
#get coverage of variable edges
with open(args.ifn_R1) as R1:
    fastq_line_count = 0
    for line in R1:
        if (fastq_line_count % 400000) == 0:
            print("Fastq 1 reads count",fastq_line_count/4)
        if fastq_line_count % 4 == 1:
            full_seq = line.rstrip()
            for base in range(0, len(full_seq)-(k_length)+1):
                k_seq = full_seq[base:base+k_length]
                k_comp_sequence = reverse_complement(k_seq) 
                if k_sequences_to_label.get(k_seq) != None:
                    k_label = k_sequences_to_label[k_seq]
                    k_coverage[k_label] += 1
                elif k_sequences_to_label.get(k_comp_sequence) != None:
                    k_label = k_sequences_to_label[k_comp_sequence]
                    k_coverage[k_label] += 1   

        fastq_line_count += 1
        
print("Now going through R2...")
#get coverage of variable edges
count = 0
with open(args.ifn_R2) as R2:
    fastq_line_count = 0
    for line in R2:
        if (fastq_line_count % 400000) == 0:
            print("Fastq 2 reads count",fastq_line_count/4)
        if fastq_line_count % 4 == 1:
            full_seq = line.rstrip()
            for base in range(0, len(full_seq)-(k_length)+1):
                count += 1
                k_seq = full_seq[base:base+k_length]
                k_comp_sequence = reverse_complement(k_seq) 
                if k_sequences_to_label.get(k_seq) != None:
                    k_label = k_sequences_to_label[k_seq]
                    k_coverage[k_label] += 1
                elif k_sequences_to_label.get(k_comp_sequence) != None:
                    k_label = k_sequences_to_label[k_comp_sequence]
                    k_coverage[k_label] += 1   
        fastq_line_count += 1
print("k count",count)
print("Total paired reads", fastq_line_count/4)

nodes = []
omatrix = open(args.ofn_edge_cov_matrix,"w+")
with open(args.ifn_adj_matrix) as adj_matrix:
    omatrix.write("\t")
    for line in adj_matrix:
        line = line.rstrip().split()
        for node in line:
            nodes.append(node)
            omatrix.write(node + "\t")
        omatrix.write("\n")
        break

print("Now writing coverages...")

ocov = open(args.ofn_edge_cov_list,"w+")
ocov.write("parent\tchild\tcoverage\tlength\tweight\n")
for parent in nodes:
    omatrix.write(parent + "\t")
    for child in nodes: 
        parent_dir = parent[-1]
        child_dir = child[-1]
        nodes_connected = connected.get((parent,child), False)
        if nodes_connected:
            if parent_dir == "-":
                parent_contig = parent[:-1]
                child_contig = child[:-1]
                contig_coverage = str(contig_to_coverage[parent_contig[:-1]])
                omatrix.write(contig_coverage + "\t")
                contig_length = str(contig_to_length[parent_contig[:2]])
                contig_weight = str(round(float(contig_length) / float(contig_coverage),3))
                ocov.write(parent + "\t" + child + "\t" + contig_coverage + "\t" + contig_length + "\t" + contig_weight + "\n") 
            else:
                edge_k = edge_to_k.get((parent,child))
                if edge_k != None:
                    edge_coverage = str(k_coverage[edge_k])
                    edge_length = str(k_length)
                    if edge_coverage == "0":
                        edge_coverage = str(0.00000001)
                    weight = str(round(float(edge_length) / float(edge_coverage),10))
                    omatrix.write(edge_coverage + "\t")
                    ocov.write(parent + "\t" + child + "\t" + edge_coverage + "\t" + str(k_length) + "\t" + weight + "\n")
        else:
            omatrix.write("0\t")
    omatrix.write("\n")