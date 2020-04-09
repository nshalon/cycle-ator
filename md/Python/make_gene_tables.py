import argparse
import numpy as np

#define arguments
parser = argparse.ArgumentParser(description='Map reads to plasmids and classify')
parser.add_argument('ifn_prodigal', metavar='<prodigal ifn>', type=str, help='Input prodigal table')
parser.add_argument('ifn_uniref', metavar='<uniref input>', type=str, help='Input uniref table')
parser.add_argument('ifn_summary', metavar='<summary table input>', type=str, help='Input summary table')
parser.add_argument('out_dir', metavar='<out_dir>', type=str, help='Output tables')

#parse arguments
args = parser.parse_args()

print("Input prodigal:",args.ifn_prodigal)
print("Input uniref:",args.ifn_uniref)
print("Output directory:",args.out_dir)

with open(args.ifn_prodigal) as prodigal:
    plas_len_to_genes = {}
    gene_to_first_coord = {}
    gene_to_last_coord = {}
    next(prodigal)   
    for line in prodigal:
        line = line.split()
        strand = line[4]
        plas_len = line[1][6:].split('_')     
        plas_len = plas_len[2]
        gene = line[0]
        plas_len_to_genes.setdefault(plas_len,[]).append(gene)
        if(strand == '+'):
            first_coord = line[2]
            last_coord = line[3]
        elif(strand == '-'):
            first_coord = line[3]
            last_coord = line[2]
        gene_to_first_coord[gene] = first_coord
        gene_to_last_coord[gene] = last_coord

with open(args.ifn_summary) as summary:
    plas_lengths = []
    next(summary)
    for line in summary:
        plas_len = line.split()[2]
        plas_lengths.append(plas_len)

ofile = open(args.out_dir,"w+")
ofile.write('PCE_num,PCE_len,gene,uniref,identity,coverage,evalue,bitscore,prot_desc,tax,uniref_count,start_coord,stop_coord\n')
plas_count = 0
for plas_len in plas_lengths:
    plas_count += 1
    genes = plas_len_to_genes.get(plas_len,'None')
    with open(args.ifn_uniref) as uniref:
        for line in uniref:
            test_gene = line.split()[0]
            line = line.rstrip('\n')
            line = line.replace(',',' ')
            line = line.replace('\t',',')
            if test_gene in genes:            
                start_coord = gene_to_first_coord.get(test_gene)
                end_coord = gene_to_last_coord.get(test_gene)
                ofile.write(str(plas_count) + ',' + plas_len + ',' + line + ',' + start_coord + ',' + end_coord + '\n')