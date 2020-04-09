import argparse
from collections import defaultdict
import sys
import heapq

#############################################################################
# Purpose: find the number of cycles whose (in+out)<bottleneck
#############################################################################

parser = argparse.ArgumentParser(description='rename reads to have FW and RV')
parser.add_argument('ifn_cycles', metavar='<summary of all cycles>', type=str, help='input cycles summary file')
parser.add_argument('ofn_good_cycles', metavar='<just output cycles that pass criteria>', type=str, help='input cycles summary file')

args = parser.parse_args()

header_keys = {}
with open(args.ifn_cycles) as cycles:
    header = cycles.readline()
    header = header.rstrip().split()
    for index,header_part in enumerate(header):
        header_keys[header_part] = index

ofile = open(args.ofn_good_cycles,'w+')
with open(args.ifn_cycles) as cycles:
    ofile.write(cycles.readline())
    next(cycles)
    cyc_count = 0
    cyc_count_nonmin = 0
    big_count = 0
    for line in cycles:
        big = False
        line_arr = line.rstrip().split()
        bn = float(line_arr[header_keys['bottleneck_cov']])
        out_c = float(line_arr[header_keys['total_out_cov']])
        in_c = float(line_arr[header_keys['total_in_cov']])
        if float(line_arr[header_keys['edge_count']]) == 4:
            big = True
        if bn > (out_c + in_c):
            ofile.write(line)
            cyc_count += 1
            if bn > 0.001:
                cyc_count_nonmin += 1
            if big:
                big_count+=1
                print(line_arr[header_keys['cycle']])

print("Number of cycles that pass criteria:", cyc_count)
print("Number of good cycles with greater than min cov:", cyc_count_nonmin)
print("Number of good cycles that are two+ contigs:", big_count)
             
            
            
    
    
            
            
            
            
        
    
        
