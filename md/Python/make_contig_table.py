import argparse

parser = argparse.ArgumentParser(description='make contig table from fasta and clean megahit fasta headers')
parser.add_argument('in_fasta', metavar='<paired_table>', type=str, help='paired table of reads')
parser.add_argument('out_table', metavar='<table with path ifn>', type=str, help='Input table w/ plasmid contigs')
parser.add_argument('out_fasta', metavar='<table with path ifn>', type=str, help='Input table w/ plasmid contigs')

args = parser.parse_args()

ofasta = open(args.out_fasta,'w+')
ofile = open(args.out_table,'w+')
ofile.write('contig\tlength\n')
with open(args.in_fasta) as fasta:
    for line in fasta:
        line = line.split()
        if line[0][0] == '>':
            ofasta.write(line[0] + '\n')
            ofile.write(line[0][1:] + '\t' + line[3][4:] + '\n')
        else:
            ofasta.write(line[0] + '\n')
