# #!/usr/bin/env python

import argparse

#####################################################################################
# command line arguments
#####################################################################################

# define arguments
parser = argparse.ArgumentParser(description='Create new megahit fasta with updated titles')
parser.add_argument('ifn_fasta', metavar='<fasta ifn>', type=str, help='Input assembly fasta')
parser.add_argument('ofn_fasta', metavar='<ofn fasta>', type=str, help='Output fasta file')
parser.add_argument('k', metavar='<k>', type=str, help='kmer size used in assembly')

args = parser.parse_args()
print ("ifn fasta: ", args.ifn_fasta)
print ("ofn fasta: ", args.ofn_fasta)

k_leng = len(args.k)

ofile = open(args.ofn_fasta,'w+')
with open(args.ifn_fasta) as ifasta:
    for line in ifasta:
        if line[0] == '>':
            line = line.split()
            contig_num = line[0][(3+k_leng):]
            contig_cov = line[2][6:]
            contig_len = line[3][4:]
            ofile.write(">EDGE_"+contig_num+"_length_"+contig_len+"_cov_"+contig_cov+"\n")
        else:
            ofile.write(line)
            