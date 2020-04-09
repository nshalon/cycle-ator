import argparse
import numpy as np
import os
import re

#define arguments
parser = argparse.ArgumentParser(description='make gene table and make the matrix')
parser.add_argument('ifn_prodigal', metavar='<prodigal ifn>', type=str, help='Input prodigal table')
parser.add_argument('in_pfam', metavar='<pfram output>', type=str, help='Input uniref table')
parser.add_argument('ifn_summary', metavar='<summary table input>', type=str, help='Input summary table')
parser.add_argument('out_table', metavar='<out_dir>', type=str, help='Output tables')
parser.add_argument('out_matrix', metavar='<out_dir>', type=str, help='Output tables')
#parse arguments
args = parser.parse_args()

print("Input prodigal:",args.ifn_prodigal)
print("Output table:",args.out_table)
print("Output matrix:",args.out_matrix)


with open(args.ifn_prodigal) as prodigal:
    plas_genes_to_len = {}
    gene_to_first_coord = {}
    gene_to_last_coord = {}
    next(prodigal)   
    for line in prodigal:
        line = line.split()
        plas_len = line[1][6:].split('_')     
        plas_len = plas_len[2]
        gene = line[0]
        plas_genes_to_len.setdefault(gene,[]).append(plas_len)
        first_coord = line[2]
        last_coord = line[3]
        gene_to_first_coord[gene] = first_coord
        gene_to_last_coord[gene] = last_coord

plas_length_to_plas = {}
with open(args.ifn_summary) as summary:
    plas_lengths = []
    next(summary)
    for line in summary:
        plas_len = line.split()[2]
        plas = line.split()[0]
        plas_lengths.append(plas_len)
        plas_length_to_plas[plas_len] = plas
        
ogenetable = open(args.out_table,'w+')
plas_to_classification = {}     
with open(args.in_pfam) as pfam:
    for line in pfam:
        if line[0] != '#':
            line = line.split()
            gene = line[2]
            name = line[0]
            leng = plas_genes_to_len[gene][0]
            plasmid = plas_length_to_plas[leng]
            plas_to_classification[plasmid] = 'AMR'
            ogenetable.write(plasmid + ',')
            for line_part in line:
                ogenetable.write(line_part + ',')
            ogenetable.write('\n')

print(plas_to_classification.keys())
omatrix = open(args.out_matrix,'w+')   
for length in plas_lengths:
    plas = plas_length_to_plas[length]
    print(plas)
    if plas in plas_to_classification.keys():
        omatrix.write(plas + '\n')   