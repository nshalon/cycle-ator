import argparse

parser = argparse.ArgumentParser(description='remove reads that map to contigs that arent part of cycles')
parser.add_argument('in_fasta', metavar='<in directory>', type=str, help='dir')
parser.add_argument('out_table', metavar='<removed cycles>', type=str, help='removed cycles')
parser.add_argument('out_fasta', metavar='<>', type=str, help='')
args = parser.parse_args()

ofile = open(args.out_table,'w+')
ofile.write('contig\tlength\tbin\tcontig.org\tstart.org\tend.org\n')
ofasta = open(args.out_fasta,'w+')
with open(args.in_fasta) as fasta:
    bin_count = 1
    for line in fasta:
        if line[0] == '>':
            split_line = line.split('_')
            length = split_line[3]
            contig_name = 'cyc_'+split_line[1]
            contig_name_main = contig_name + ':s1'
            ofasta.write('>'+contig_name_main+'\n')
            ofile.write(contig_name_main+'\t'+length+'\t'+str(bin_count)+'\t'+contig_name+'\t'+str(1)+'\t'+str(length)+'\n')
            bin_count += 1
        else:
            ofasta.write(line)
