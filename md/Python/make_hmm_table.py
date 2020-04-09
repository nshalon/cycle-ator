import argparse
import numpy as np
import os

#define arguments
parser = argparse.ArgumentParser(description='make gene table and make the matrix')
parser.add_argument('ifn_prodigal', metavar='<prodigal ifn>', type=str, help='Input prodigal table')
parser.add_argument('in_dir', metavar='<hmm tables>', type=str, help='Input uniref table')
parser.add_argument('ifn_summary', metavar='<summary table input>', type=str, help='Input summary table')
parser.add_argument('out_table', metavar='<out_dir>', type=str, help='Output tables')
parser.add_argument('out_matrix', metavar='<out_dir>', type=str, help='Output tables')

#parse arguments
args = parser.parse_args()

print("Input prodigal:",args.ifn_prodigal)
print("Output table:",args.out_table)
print("Output matrix:",args.out_matrix)

with open(args.ifn_prodigal) as prodigal:
    plas_len_to_genes = {}
    gene_to_first_coord = {}
    gene_to_last_coord = {}
    next(prodigal)   
    for line in prodigal:
        line = line.split()
        plas_len = line[1][6:].split('_')     
        plas_len = plas_len[2]
        gene = line[0]
        plas_len_to_genes.setdefault(plas_len,[]).append(gene)
        first_coord = line[2]
        last_coord = line[3]
        gene_to_first_coord[gene] = first_coord
        gene_to_last_coord[gene] = last_coord

dirs = os.listdir(args.in_dir)
gene_to_profile = {}
for hmm in dirs:
    hmm_dir = args.in_dir + '/' + hmm
    with open(hmm_dir) as hmm_in:
        for line in hmm_in:
            if line[0] == '#':
                continue
            else:
                line = line.split()
                gene = line[2]
                protein = line[0]
                if '-' not in line[4]:
                    eval = '2'
                else:
                    eval = line[4].split('-')[1]
                protein_eval = protein + ',' + eval
                gene_to_profile.setdefault(gene,[]).append(protein_eval)

plas_length_to_plas = {}
with open(args.ifn_summary) as summary:
    plas_lengths = []
    next(summary)
    for line in summary:
        plas_len = line.split()[2]
        plas = line.split()[0]
        plas_lengths.append(plas_len)
        plas_length_to_plas[plas_len] = plas
        
otable = open(args.out_table,'w+')
otable.write("cycle,cycle_len,gene,gene_start,gene_stop,hmm_protein,eval\n")
for plas_length in plas_lengths:
    plas = plas_length_to_plas[plas_length]
    genes = plas_len_to_genes.get(plas_length,'None')
    if genes == 'None':
        continue
    for gene in genes:
        gene_start = gene_to_first_coord[gene]
        gene_stop = gene_to_last_coord[gene]
        profiles = gene_to_profile.get(gene,'None')
        if profiles == 'None':
            continue
        for profile in profiles:
            protein = profile.split(',')[0]
            eval_exp = profile.split(',')[1]
            eval = '1*10^-' + eval_exp
            otable.write(plas + ',' + plas_length + ',' + gene + ',' + gene_start + ',' + gene_stop + ',' + protein + ',' + eval + '\n') 

protein_to_index = {}
protein_ordered = []
index = 0
for dir in dirs:
    index += 1
    protein = dir[:-4]
    protein_to_index[protein] = index
    protein_ordered.append(protein)
    
omatrix = open(args.out_matrix,'w+')
for protein in protein_ordered:
    omatrix.write(protein+'\t')
omatrix.write('\n')


for plas_length in plas_lengths:
    plas = plas_length_to_plas[plas_length]
    genes = plas_len_to_genes.get(plas_length,'None')
    protein_evals = [0] * len(dirs)
    if genes == 'None':
        for eval in protein_evals:
            omatrix.write(str(eval) + '\t')
        omatrix.write('\n')
        continue
    for gene in genes:
        profiles = gene_to_profile.get(gene,'None')
        if profiles == 'None':
            continue
        for profile in profiles:
            protein = profile.split(',')[0]
            eval_exp = profile.split(',')[1]
            protein_evals[protein_to_index[protein]] = eval_exp            
    for eval in protein_evals:
        omatrix.write(str(eval) + '\t')
    omatrix.write('\n')

        
        
    
    
                
    
        
