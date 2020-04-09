#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='Create buckets for classification table')
parser.add_argument('ifn_table', metavar='<table>', type=str, help='Input classified table')
parser.add_argument('binsize', metavar='<size of bin>', type=int, help='Input desired size of bin')
parser.add_argument('ofn_bintable', metavar='<ofn>', type=str, help='Output dir. for binned table')
# parse arguments
args = parser.parse_args()
print ("ifn classified table: ", args.ifn_table)
print ("binsize: ", args.binsize)
print ("ofn bintable: ", args.ofn_bintable)

print("Writing classified table with bins the size of:",args.binsize)

endcount = 0
for line in open(args.ifn_table).readlines(  ): endcount += 1
endcount -= 1
ofile = open(args.ofn_bintable, "w+")
with open(args.ifn_table) as ifile:
    linecount = 0
    next(ifile)
    ofile.write("Plasmid\tBase_Pair\tStrand\tGood_Reads\tDistant_Reads\tNot_on_plasmid\tOn_lin_seq\n")
    bin = args.binsize
    count = 0
    one = two = three = four = 0
    strand = None
    for line in ifile:
        linecount+=1
        line = line.split()
        if ((strand != line[2] or linecount==endcount) and count != 0):
            scale = bin/count
            if(linecount==endcount):
                ofile.write(plasmid + '\t' + str(int(bp)+1) + '\t' + strand + '\t' + str(one*scale) + '\t' + str(two*scale) + '\t' + str(three*scale) + '\t' + str(four*scale) + '\n')
            else:
                ofile.write(plasmid+'\t'+bp+'\t'+strand+'\t'+str(scale*one)+'\t'+str(scale*two)+'\t'+str(scale*three)+'\t'+str(scale*four)+'\n')
            one = two = three = count = 0
        plasmid = line[0]
        bp = line[1]
        strand = line[2]
        one += int(line[3])
        two += int(line[4])
        three += int(line[5])
        four += int(line[6])
        count+=1
        if count >= bin:
            ofile.write(plasmid+'\t'+bp+'\t'+strand+'\t'+str(one)+'\t'+str(two)+'\t'+str(three)+'\t'+str(four)+'\n')
            one = two = three = count = four = 0
print("Find results @",args.ofn_bintable)
