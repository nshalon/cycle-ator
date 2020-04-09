import argparse
import numpy as np


parser = argparse.ArgumentParser(description='Calculate read length distribtion and calculate coverage')
parser.add_argument('ifn_pairedtable', metavar='<paired_table>', type=str, help='paired table of reads')
parser.add_argument('ifn_summary', metavar='<table with path ifn>', type=str, help='Input table w/ plasmid contigs')
parser.add_argument('ofn_covstats', metavar = '<ofn_covstats>', type=str, help='Output coverage and read length stats')

args = parser.parse_args()

with open(args.ifn_summary) as summary:
    next(summary)
    plas_tables = []
    for line in summary:
        line = line.split()
        length = int(line[2])
        plas_table = [0] * length
        plas_tables.append(plas_table)

col_htable = {}
with open(args.ifn_pairedtable) as itable:
    line = itable.readline()
    header = line.split()
    index = 0
    for col in header:
        col_htable[col] = index
        index += 1

with open(args.ifn_pairedtable) as table:
    next(table)
    for line in table:
        line = line.split()
        plas = line[col_htable['contig']].split(':')[0]
        plas_num = int(plas[4:]) - 1
        coords = []
        coords.append(int(line[col_htable['coord']]))
        coords.append(int(line[col_htable['back_coord']]))
        first_coord = min(coords[0],coords[1])
        last_coord = max(coords[0],coords[1])
        for i in range(first_coord-1,last_coord):
            plas_tables[plas_num][i] += 1

ofile = open(args.ofn_covstats,'w+')
ofile.write('cycle\tcoord\tcount\n')
table_count = 0
for plas_num in range(0,len(plas_tables)):
    table = plas_tables[plas_num]
    length = len(table)
    table_count += 1
    for i in range(0,length):
        if (str(table[i])!= '0'):
            print(str(table[i]))
        ofile.write(str(table_count) + '\t' + str(i) + '\t' + str(table[i]) + '\n')
        
            
        