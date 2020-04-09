import argparse
from collections import defaultdict
import sys
import heapq

#############################################################################
# Purpose: find the number of cycles whose (in+out)<bottleneck
#############################################################################

parser = argparse.ArgumentParser(description='rename reads to have FW and RV')
parser.add_argument('ifn_reads', metavar='<summary of all cycles>', type=str, help='input cycles summary file')

args = parser.parse_args()

with open(args.ifn_reads) as reads:
    readc = -1
    ncount = 0
    small_n = 0
    sense_reads = 0
    for line in reads:
        if readc % 4 == 0:
            seq = line.rstrip()
            if "N" in seq:
                ncount += 1
            nseq = [1 for c in seq if c == "N"]
            num_n = len(nseq)
            if num_n < 7:
                small_n += 1
            if num_n < 50:
                sense_reads += 1
        readc += 1

num_reads = readc/4
print("Total reads:",num_reads)
print("Total reads with N:",ncount)
print("Percent reads with N:", float(ncount)/num_reads ) 
print("Percent reads that have usable info:",sense_reads/num_reads)
print("Percent reads that have less than 7 N's:",small_n/num_reads)
print("Percent of usable reads with less than 7 N's:",small_n/sense_reads)


