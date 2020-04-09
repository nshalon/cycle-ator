#!/usr/bin/env python
import argparse
import re

parser = argparse.ArgumentParser(description='Create new fasta w/ renamed contigs')
parser.add_argument('ifn_contigfasta', metavar='<contigfasta>', type=str, help='the renamed contig fasta')
parser.add_argument('ifn_plasfasta', metavar='<plasmidfasta>', type=str, help='the fasta of all plasmids')
parser.add_argument('ifn_pairedtable', metavar='<paired_table>', type=str, help='paired table of reads')
parser.add_argument('i_path', metavar='<input_table_path>', type=str, help='input table path')
parser.add_argument('o_path', metavar='<output_table_path>', type=str, help='output table path')
parser.add_argument('k', metavar='<k>', type=int, help='kmer size used in assembly')
# parse arguments
args = parser.parse_args()

print ("ifn contigfasta: ", args.ifn_contigfasta)
print ("ifn plasfasta: ", args.ifn_plasfasta)
print ("ifn paired tabled", args.ifn_pairedtable)
print ("ifn path: ", args.i_path)
print ("ofn path: ", args.o_path)

print("Finding orientations of contigs within plasmid")

k = args.k

ofile = open(args.o_path, "w+")
ofile.write("Plasmid\tIndex\tContig\tContig_len\tCumsum\tOrientation\n")
with open(args.i_path, 'r') as path:
    plasmid = None
    next(path)
    search = None
    contiglist = [None,None]
    for line in path:
        ofile.write(line.rstrip())
        line = line.split()
        contiglist[0]=contiglist[1]
        contiglist[1]=line[2]
        print(contiglist)
        if plasmid != line[0]:
            plasmid = line[0]
            contig = line[2]
            with open(args.ifn_contigfasta) as contigs:
                print(contig)
                searchline = False
                for line in contigs:
                    if(searchline):
                        if len(line) < (k+20):
                            search = line[k:]
                        else:
                            search = line[k:k+20]
                            print('There was search')
                        break
                    if contig in line:
                        searchline = True
            with open(args.ifn_plasfasta) as plasmids:
                plasmidnum = int(plasmid[1:])
                searchfasta = False
                for line in plasmids:
                    if(searchfasta): 
                        if search in line[:100]:
                            orientation = '+'
                        else:
                            orientation = '-'
                        break
                    if line[0] == '>':
                        headernum=re.sub('_length.+', "", line)
                        headernum=int(headernum[7:])
                    if headernum == plasmidnum:
                        searchfasta = True

        else:
            table = open(args.ifn_pairedtable)
            for line in table:
                line = line.split()
                if line[1]==contiglist[0] and line[14] == contiglist[1]:
                    if line[4] == line[17]:
                        if orientation == '+':
                            orientation = '-'
                        else:
                            orientation = '+'
                    break
        ofile.write('\t'+ orientation + '\n')