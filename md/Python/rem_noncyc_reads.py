import argparse

#define arguments
parser = argparse.ArgumentParser(description='remove reads that map to contigs that arent part of cycles')
parser.add_argument('ifn_pairedtable', metavar='<ifn_paired table reads', type=str, help='Input paired table of reads')
parser.add_argument('ofn_pairedtable', metavar='<ifn_paired table reads w/ reads removed>', type=str, help='Output table')
parser.add_argument('ofn_stats', metavar='<stats>', type=str, help='stats')

#parse arguments
args = parser.parse_args()

print('Now taking out irrelevant reads...')
otable = open(args.ofn_pairedtable,'w+')
line_count = 0
col_htable = {}
write_count = 0
index = 0
with open(args.ifn_pairedtable) as itable:
    for read in itable:
        if(line_count == 0):
            header = read.split()
            for col in header:
                col_htable[col] = index
                index += 1
            otable.write(read)
            line_count += 1
            continue
        line_count += 1
        line = read.split()
        if line[col_htable['contig1']][0] != 'C' or line[col_htable['contig2']][0] != 'C':
            otable.write(read)
            write_count += 1
            
ofile = open(args.ofn_stats,'w+')
ofile.write('Write count %s \n' % write_count)
ofile.write('Total count %s \n' % line_count)
print('Write count %s' % write_count)
print('Total count %s' % line_count)
print('Reads in cycle is %s percent of total reads' % (str(write_count/line_count)))