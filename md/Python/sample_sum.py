#!/usr/bin/env python

import argparse
import re
import numpy as np

parser = argparse.ArgumentParser(description='Create new fasta w/ renamed contigs')
parser.add_argument('ifn_CV', metavar='<input CV table>', type=str, help='the table with all plasmid CV vals')
parser.add_argument('ifn_plasmidpath', metavar='<plasmid path table>', type=str, help='the table of contigs in plasmid')
parser.add_argument('ifn_bin_summary', metavar='<binned summary table>', type=str, help='summary table after binning with 100')
parser.add_argument('ofn_summary', metavar='<output summary table>', type=str, help='output summary table')
parser.add_argument('k', metavar='<k>', type=int, help='kmer size used in assembly')
# parse arguments
args = parser.parse_args()

print ("ifn CV table: ", args.ifn_CV)
print ("ifn plasmidpath: ", args.ifn_plasmidpath)
print ("ifn bin summary", args.ifn_bin_summary)
print ("ofn summary path: ", args.ofn_summary)

k = args.k

ofile = open(args.ofn_summary,"w+")
ofile.write('Plasmid\tNumber_of_contigs\tLength\tCV\tMed_density\tMin_density\n')
with open(args.ifn_CV) as CV:
    next(CV)
    for line in CV:
        line = line.split()
        plasmid = line[0]
        CV_val = round(float(line[1]),3)
        with open(args.ifn_bin_summary) as bin:
            append = False
            supporting_reads = []
            for binline in bin:
                binline = binline.split()
                if binline[0] == plasmid:
                    append = True
                if append:
                    if binline[0] != plasmid:
                        break
                    supporting_reads.append(int(float((binline[3]))))
            med_density = round(np.median(supporting_reads)/100,3)
            if(len(supporting_reads)!=0):
                min_density = round(np.min(supporting_reads)/100,3)
        index_num = 0
        with open(args.ifn_plasmidpath) as path:
            length = 0
            right_plasmid = False
            for pathline in path:
                pathline = pathline.split()
                if right_plasmid and pathline[0] != plasmid_in_path:
                    break
                if pathline[0][1:] == plasmid:
                    right_plasmid = True
                    plasmid_in_path = pathline[0]
                    index_num = pathline[1]
                    length = int(pathline[4]) - k + int(pathline[3])
        
        ofile.write(plasmid + '\t' + str(index_num) + '\t' + str(length) + '\t' + repr(CV_val) + '\t' + repr(med_density) + '\t' + repr(min_density) + '\n')
