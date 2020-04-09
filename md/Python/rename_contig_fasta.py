import argparse

#############################################################################
# Purpose: input a fasta (of contigs) and file showing how each original
# contig name is mapped to the renamed name. Then, output a fasta with renamed
# headers
#############################################################################

parser = argparse.ArgumentParser(description='rename reads to have FW and RV')
parser.add_argument('ifn_fasta', metavar='<contig fasta input to be renamed>', type=str, help='assembly output of contigs')
parser.add_argument('ifn_name_map', metavar='<file with original names in first column, renamed name in second>', type=str, help='Original_contig new_contig')
parser.add_argument('ofn_rename_fasta', metavar='<renamed fasta>', type=str, help='renamed fasta')

args = parser.parse_args()

old_contig_to_new = {}
with open(args.ifn_name_map) as map:
    for line in map:
        line = line.split()
        contig_num = line[0].split("_")[1]
        old_contig_to_new[contig_num] = line[1]
        
ofile = open(args.ofn_rename_fasta,'w+')
with open(args.ifn_fasta) as fasta:
    line_count = 0
    keyError = False
    for line in fasta:
        if keyError:
            keyError = False
            continue
        if line[0] == ">":
            if line_count % 1000 == 0:
                print("contig",line_count)
            line_count += 1
            line = line.split()
            
            contig_num = line[0].split("_")[1]
            if old_contig_to_new.get(contig_num) == None:
                keyError = True
                continue
                
            new_name = old_contig_to_new[contig_num]
            ofile.write(">" + new_name + "\n")
        else:
            ofile.write(line)
            
args = parser.parse_args()

old_contig_to_new = {}
with open(args.ifn_name_map) as map:
    for line in map:
        line = line.split()
        contig_num = line[0].split("_")[1]
        old_contig_to_new[contig_num] = line[1]
        
ofile = open(args.ofn_rename_fasta,'w+')
with open(args.ifn_fasta) as fasta:
    line_count = 0
    keyErrorCount = 0
    keyError = False
    for line in fasta:
        if keyError:
            keyError = False
            continue
        if line[0] == ">":
            if line_count % 1000 == 0:
                print("contig",line_count)
            line_count += 1
            line = line.split()
            
            contig_num = line[0].split("_")[1]
            if old_contig_to_new.get(contig_num) == None:
                keyError = True
                keyErrorCount += 1
                continue
                
            new_name = old_contig_to_new[contig_num]
            ofile.write(">" + new_name + "\n")
        else:
            ofile.write(line)
            
print("Number of key errors:",keyErrorCount)
            