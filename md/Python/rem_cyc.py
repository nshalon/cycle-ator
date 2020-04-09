import argparse
import os
import re
import numpy as np

parser = argparse.ArgumentParser(description='remove reads that map to contigs that arent part of cycles')
parser.add_argument('in_dir', metavar='<in directory>', type=str, help='dir')
parser.add_argument('rem_cyc', metavar='<removed cycles>', type=str, help='removed cycles')
args = parser.parse_args()

dirs = os.listdir(args.in_dir)

print('In dir:',args.in_dir)

subject_dirs = []
for dir in dirs:
    if re.match('s[0-9]+_',dir):
        subject_dirs.append(dir)

print(subject_dirs)

cyc_table = []
for dir in subject_dirs:
    cyc_count = 0
    file = args.in_dir + '/' + dir + '/stats'
    with open(file) as stats:
        for line in stats:
            cyc_count += 1
    print(cyc_count)
    cyc_table = [0] * (cyc_count-1)
    file = args.in_dir + '/' + dir + '/cov_table'
    file = open(file)
    line = file.readline()
    header = line.split()
    line = line.rstrip()
    index = 0
    col_htable = {}
    for col in header:
        col_htable[col] = index
        index += 1
    break

for dir in subject_dirs:
    file = args.in_dir + '/' + dir + '/cov_table'
    with open(file) as cov_table:
        cyc = '1'
        good_cyc = True
        next(cov_table)
        for line in cov_table:
            line = line.split()
            if(cyc != line[col_htable['Plasmid']]): #fix so last cyc gets tested
                if good_cyc:
                    cyc_table[int(cyc)-1] += 1
                good_cyc = True
            cyc = line[col_htable['Plasmid']]
            cov = int(line[col_htable['Cov']])
            out_cov_pos = int(line[col_htable['Out_cyc_pos']])
            out_cov_neg = int(line[col_htable['Out_cyc_neg']])
            if cov < 2 or cov<out_cov_pos or cov<out_cov_neg:
                good_cyc = False

ofile = open(args.rem_cyc,'w+')
for cyc in range(0,len(cyc_table)):
    if cyc_table[cyc] > 5:
        print('good',str(cyc+1))
    if cyc_table[cyc] == 0:
        print('removed:',str(cyc+1))
        ofile.write(str(cyc+1)+'\n')        


        