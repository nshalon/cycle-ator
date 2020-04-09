import argparse
import os
import re
import numpy as np
import math

parser = argparse.ArgumentParser(description='remove reads that map to contigs that arent part of cycles')
parser.add_argument('in_dir', metavar='<in directory>', type=str, help='dir')
parser.add_argument('med_cov_table', metavar='<removed cycles>', type=str, help='removed cycles')
args = parser.parse_args()

dirs = os.listdir(args.in_dir)

print('In dir:',args.in_dir)

subject_dirs = []
for dir in dirs:
    if re.match('s[0-9]+_',dir):
        subject_dirs.append(dir)
        
subj_dirs_ordered = []
for i in range(0,2*len(subject_dirs)):
    print('hello')
    for j in range(1,3):
        for dirs in subject_dirs:
            dirs_test = dirs.split('_')
            if int(dirs_test[0][1:]) == i and int(dirs_test[1][1:]) == j:
                subj_dirs_ordered.append(dirs)
                print('hello')
for i in range(0,len(subject_dirs)):
    subject_dirs[i] = subj_dirs_ordered[i]
    print(subj_dirs_ordered[i])

for dir in subject_dirs:
    cyc_count = 0
    file = args.in_dir + '/' + dir + '/stats'
    with open(file) as stats:
        for line in stats:
            if cyc_count == 0:
                index = 0
                col_htable = {}
                line = line.split()
                for col in line:
                    col_htable[col] = index
                    index += 1
            cyc_count += 1
    cov_table = np.zeros((len(subject_dirs), cyc_count-1), dtype=int)
    break


subject_num=0
for dir in subject_dirs:
    file = args.in_dir + '/' + dir + '/stats' 
    with open(file) as stats:
        next(stats)
        cyc_count = 0
        for line in stats: 
            line = line.split()
            med_cov = float(line[col_htable['Med_cov']])
            med_cov = math.floor(med_cov)
            cov_table[subject_num][cyc_count] = med_cov
            cyc_count += 1
    subject_num += 1

ofile = open(args.med_cov_table,'w+')
for sub in subject_dirs:
    ofile.write(sub + '\t')
ofile.write('\n')


# for subj_count in range(0,len(subject_dirs)):
#     print(cov_table[subj_count])

for cyc_count in range(0,len(cov_table[0])):  
    for subj_count in range(0,len(subject_dirs)):
        cyc_counter = 0
        for cov in cov_table[subj_count]:
            if cyc_counter == cyc_count:
                ofile.write(str(cov) + '\t')
            cyc_counter += 1
    ofile.write('\n')
    
    
    