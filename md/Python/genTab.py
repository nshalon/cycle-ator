import argparse

# define arguments below
# each run takes in a filename of all the original reads (preferable in .gz format)
# and outputs a list of unique name (minus FW or RV .fastq.gz) and each read name (2 read names, one FW and one RV, per unique name)

parser = argparse.ArgumentParser(description='rename reads and output edited unique names and all edited read names')
parser.add_argument('ifn_fa', metavar='<fasta file to create a tab file>', type=str, help='')
parser.add_argument('ofn_tab', metavar='<tab file>', type=str, help='tab file')

args = parser.parse_args()

ofile = open(args.ofn_tab,'w+')
ofile.write("contig\tlength\txcov\tcirc\n")
with open(args.ifn_fa) as fasta:
    size = 0
    linenum = 0
    name = ''
    for line in fasta:
        if linenum == 0:
            name = line[1:].rstrip()
        linenum += 1
        size += len(line.rstrip())
    print(name,"has size",size)
    ofile.write("m1\t" + str(size) + "\t50\tTRUE\n")
        
        
