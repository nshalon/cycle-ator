import argparse
import numpy as np
from operator import itemgetter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#############################################################################
# Purpose: input a fasta table that contains the parsed output of bwa and calculate xcoverage of each contig
# Decisions:
# 1) When both sides of a read fall on same contig then add one to the coverage for bp from first tail to second tail of read
# 2) When only side of a read fall on same contig then only add coverage to the bp that fall directly on that one read side (from tail to head)
# 3) When done with calculating coverage along all contigs, take the median base coverage of that contig and divide it by 
# the length of that contig to get the xcoverage stat
# 
#############################################################################

parser = argparse.ArgumentParser(description='rename reads to have FW and RV')
parser.add_argument('ifn_bwa_table', metavar='<contig fasta input to be renamed>', type=str, help='assembly output of contigs')
parser.add_argument('ifn_name_map', metavar='<file with original names in first column, renamed name in second, length in third>', type=str, help='Original_contig new_contig')
parser.add_argument('ofn_contig_coverage', metavar='<renamed fasta>', type=str, help='renamed fasta')
parser.add_argument('ofn_contig_coverage_whole_molecules_dir', metavar='<renamed fasta>', type=str, help='renamed fasta')

args = parser.parse_args()


# connect contig name to length
# initialize a 0-array that is the length of a contig connected to the name of a contig
contig_to_cov_table = {}
with open(args.ifn_name_map) as contigs:
    for line in contigs:
        line = line.split()
        contig = line[1]
        length = int(line[2])
        contig_to_cov_table[contig] = np.zeros((length,), dtype=int)

map_col_index = {}
with open(args.ifn_bwa_table) as mapfile:
    line = mapfile.readline()
    header = line.split()
    index = 0
    for col in header:
        map_col_index[col] = index
        index += 1

with open(args.ifn_bwa_table) as mapfile:
    next(mapfile)
    linenum = 1
    for line in mapfile:
        linenum += 1
        pair = line.split()
        strand1 = pair[map_col_index['strand1']]
        strand2 = pair[map_col_index['strand2']]
        contig1 = pair[map_col_index['contig1']]
        contig2 = pair[map_col_index['contig2']]
        if (contig1 == contig2) & (strand1 == strand2):
            print(contig1,contig2,strand1,strand2)
            continue
        head1 = int(pair[map_col_index['coord1']])
        head2 = int(pair[map_col_index['coord2']])
        tail1 = int(pair[map_col_index['back_coord1']])
        tail2 = int(pair[map_col_index['back_coord2']])
        if contig1 == contig2:
            if strand1 == '1': 
                start_coord = tail1-1
                stop_coord = tail2
            else: 
                start_coord = tail2-1
                stop_coord = tail1
            for bp in range(start_coord,stop_coord):
                contig_to_cov_table[contig1][bp] += 1
        else:
            if strand1 == '1':
                start_coord = tail1 - 1
                stop_coord = head1
            else:
                start_coord = head1 - 1
                stop_coord = tail1
            for bp in range(start_coord,stop_coord):
                contig_to_cov_table[contig1][bp] += 1
            if strand2 == '1':
                start_coord = tail2 - 1
                stop_coord = head2
            else:
                start_coord = head2 - 1
                stop_coord = tail2
            for bp in range(start_coord,stop_coord):
                contig_to_cov_table[contig2][bp] += 1


ographdir = args.ofn_contig_coverage_whole_molecules_dir
ofile = open(args.ofn_contig_coverage,'w+')
contig_covs = [(k,v) for k,v in contig_to_cov_table.items()]
sorted_contig_covs = sorted(contig_covs, key=itemgetter(0))
ofile.write("contig coverage\n")


for item in sorted_contig_covs:
    contig = item[0]
    cov_arr = contig_to_cov_table[contig]
    contig_length = len(cov_arr)
    median_cov = np.sum(cov_arr)/contig_length
    ofile.write(contig + " " + str(median_cov) + "\n")
    hist = plt.plot(range(contig_length),cov_arr)
    plt.title(contig + "coverage distribution")
    plt.xlabel('base')
    plt.ylabel('coverage')
    ofn_graph = ographdir + '/' + contig + '_cov.pdf'
    plt.savefig(ofn_graph)
    plt.close()

                
            
                        