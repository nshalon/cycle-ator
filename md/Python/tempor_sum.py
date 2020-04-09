import argparse

#define arguments
parser = argparse.ArgumentParser(description='remove reads that map to contigs that arent part of cycles')
parser.add_argument('ifn_tablepath', metavar='<ifn_paired table reads', type=str, help='Input paired table of reads')
parser.add_argument('ofn_summ', metavar='<stats>', type=str, help='stats')
parser.add_argument('k', metavar='<k>', type=int, help='kmer size used in assembly')

#parse arguments
args = parser.parse_args()
ofile = open(args.ofn_summ,'w+')
num_lines = 0
k = args.k
for line in open(args.ifn_tablepath):
    num_lines += 1
num_lines -= 1
with open(args.ifn_tablepath) as path:
    line_count = 0
    ofile.write('Plasmid\tNum_contigs\tLength\n')
    plasmid = 'P1'
    num_contigs = 0
    plas_length = 0
    next(path)
    for line in path:
        line_count += 1
        line = line.split()
        print(num_lines,line_count,plasmid)
        if plasmid != line[0]:
            ofile.write(plasmid+'\t'+str(num_contigs)+'\t'+str(plas_length)+'\n')
            plasmid = line[0]
            plas_length = 0
            num_contigs = 0
        plasmid = line[0]
        plas_length += int(line[3]) - k
        num_contigs += 1
        if num_lines == line_count:
            ofile.write(line[0]+'\t'+str(num_contigs)+'\t'+str(plas_length)+'\n')
        