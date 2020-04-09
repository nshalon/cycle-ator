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
parser.add_argument('out_gene_table', metavar='<out_dir>', type=str, help='Output tables')
parser.add_argument('out_cycle_class', metavar='<out_dir>', type=str, help='Output tables')

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
        strand = line[4]
        plas_len = line[1][6:].split('_')     
        plas_len = plas_len[2]
        gene = line[0]
        plas_genes_to_len[gene] = plas_len
        if(strand == '+'):
            first_coord = line[2]
            last_coord = line[3]
        elif(strand == '-'):
            first_coord = line[3]
            last_coord = line[2]
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
         
raw_transposon = ["transpos.*","insertion","resolv.*","Tra[A-Z]","Tra[0-9]","IS[0-9]",]
raw_plasmid = ["resolv.*", "relax.*", "conjug.*","trb","mob.*","plasmid.*","T4SS","toxin*","partitioning","segregation","recombination","replication","pil.*","Immunity"]
raw_phage = ["capsid","phage.*","tail","head","tape","antitermination","virus.*","Bacteriophage","viral.*","sipho*","Baseplate","T4.*","myovir.*"]

reg_transposon = "(" + ")|(".join(raw_transposon) + ")"
reg_plasmid = "(" + ")|(".join(raw_plasmid) + ")"
reg_phage = "(" + ")|(".join(raw_phage) + ")"

ofile = open(args.out_table,'w+')
ogenetable = open(args.out_gene_table,'w+')
plas_to_classification = {}
ogenetable.write('PCE_num,class,start_coord,stop_coord,prot_desc,')
with open(args.in_pfam) as pfam:
    next(pfam)
    for line in pfam:
        header_parts = line.split()
        for header_part in header_parts[3:]:
            ogenetable.write(header_part + ',')
        break 
    ogenetable.write('\n')
  
with open(args.in_pfam) as pfam:
    for line in pfam:
        if line[0] != '#':
            print_line = line
            line = line.split()
            gene = line[2]
            name = line[0]
            description_parts = line[18:]
            classification = 'None'
            for description in description_parts: 
                if re.match(reg_transposon,description,flags=re.IGNORECASE):
                    classification = 'transposon'
                if re.match(reg_plasmid,description,flags=re.IGNORECASE):
                    classification = 'plasmid'
                if re.match(reg_phage,description,flags=re.IGNORECASE):
                    classification = 'phage'
            if re.match("T4SS.*",line[0]):
                classification = 'plasmid'
            length = plas_genes_to_len[gene]
            plasmid = plas_length_to_plas[length]
            plas_to_classification.setdefault(plasmid,[]).append(classification)
            gene_start = gene_to_first_coord[gene]
            gene_end = gene_to_last_coord[gene]
            if classification != 'None':
                ofile.write(plasmid + ',' + gene + ',' + name + ',' + classification + '\n')
            ogenetable.write(plasmid[1:] + ',' + classification + ',' + gene_start + ',' + gene_end + ',')
            for line_part in line:
                ogenetable.write(line_part + ',')
            ogenetable.write('\n')

omatrix = open(args.out_matrix,'w+')   
for length in plas_lengths:
    plas = plas_length_to_plas[length]
    if plas in plas_to_classification.keys():
        classifications = plas_to_classification[plas]
        omatrix.write(plas)
        for classification in classifications:
            if classification != 'None':
                omatrix.write(',' + classification)
        omatrix.write('\n')     

cycle_class = open(args.out_cycle_class,'w+')
cycle_class.write('cycle,class\n')
for length in plas_lengths:
    plas = plas_length_to_plas[length]
    cycle_class.write(plas+',')
    if plas in plas_to_classification.keys():
        classifications = plas_to_classification[plas]
    else:
        classifications = ['cryptic']
    cycle_type = 'unknown'
    if 'transposon' in classifications:
        cycle_type = 'transposon'
    if 'plasmid' in classifications:
        cycle_type = 'plasmid'
    if 'phage' in classifications:
        cycle_type = 'phage'
    if 'cryptic' in classifications:
        cycle_type = 'cryptic'
    cycle_class.write(cycle_type+'\n')
    
        
    
    
        