#!/usr/bin/env python

import argparse
import re
#####################################################################################
# command line arguments
#####################################################################################

# define arguments
parser = argparse.ArgumentParser(description='Create new fasta w/ renamed contigs')
parser.add_argument('ifn_recycler', metavar='<recycler ifn>', type=str, help='Input recycler file')
parser.add_argument('ifn_fasta', metavar='<fasta ifn>', type=str, help='Input assembly fasta')
#parser.add_argument('ifn_fastg', metavar='<fastg ifn>', type=str, help='Input assembly fastg')
parser.add_argument('ofn_fasta', metavar='<ofn fasta>', type=str, help='Output fasta file')
parser.add_argument('ofn_table', metavar='<ofn table>', type=str, help='Output table file')
parser.add_argument('k', metavar='<k>', type=int, help='kmer size used in assembly')
#parser.add_argument('ofn_fastg', metavar='<ofn fastg>', type=str, help='Output fastg file')
# parser.add_argument('plasmid', metavar='<plasmid>', type=str, help='Plasmid (for output directories for table and fasta)')
# parse arguments
args = parser.parse_args()
print ("ifn recycler: ", args.ifn_recycler)
print ("ifn fasta: ", args.ifn_fasta)
#print ("ifn fastg: ", args.ifn_fastg)
print ("ofn fasta: ", args.ofn_fasta)
print ("ofn table: ", args.ofn_table)
#print ("ofn fastg: ", args.ofn_fastg)

#####################################################################################
# add plasmid-associated contigs to hashtable
#####################################################################################
k = args.k
# hashtable with contigs as keys
htable = {}

print ('reading recycler file:',args.ifn_recycler)
with open(args.ifn_recycler) as ifile:
    ncontigs = 0
    plasmidmax = 1
    repeatnum = 1
    for line in ifile:
        line = line.rstrip("\n")
        if (line.__contains__('RNODE')):
            plasmidnum = re.sub('_length.+', "", line)
            print(plasmidnum)
            plasmidnum = int(plasmidnum[6:])
            if (plasmidmax<plasmidnum):
                plasmidmax=plasmidnum
        if (line[0]=='(' or line[0] == 'N'):
            line = line.replace("(", "").replace(",", "").replace("'", "").replace(")", "")
            line = line.replace("\"", "")
            contigs = line.split(" ")
            plasmidcontigcount = 1
            for contig in contigs:
                contig = contig[4:]
                contig = re.sub('_ID_[0-9]+$', "", contig)
                contig = re.sub('^_[0-9]+_', "", contig)
                print(contig)
                if contig in htable.keys() and htable.get(contig)[0]!= 'R':
                    htable[contig] = ('R'+str(repeatnum))
                    plasmidcontigcount += 1
                    repeatnum += 1
                elif(contig not in htable.keys()):
                    htable[contig] = ('P'+str(plasmidnum)+'_'+str(plasmidcontigcount))
                    plasmidcontigcount+=1
                if (re.match('^l',contig)):
                    ncontigs+=1


# for i in range(1,plasmidmax+1): #Puts the right suffix on the non-repeat plasmid contigs
#     plasmidcontig=1
#     for key in htable.keys():
#         if(re.match(('P'+str(i)+'_$'),htable.get(key))):
#             htable[key] = htable.get(key)+str(plasmidcontig)
#             print(htable.get(key))
#             plasmidcontig+=1

print ('number of contigs in plasmids: ',ncontigs)

# #####################################################################################
# # go over input fasta
# #####################################################################################

print(htable.keys())
print ('reading fasta file:',args.ifn_recycler)
print ('writing fasta file:',args.ofn_fasta)
ofile = open(args.ofn_fasta,"w+")
total = 0
contig = 'none'
notplasmidcount = 1
notplasmidname = {}
with open(args.ifn_fasta) as ifile:
    for line in ifile:
        line = line.rstrip("\n")
        if (line[0:1] == ">"):
            total+=1
            print('Contig: ' + line)
            contig = line[5:]
            print('Contig-pre: ' + contig)
            contig = re.sub('^_[0-9]+_', "", contig)
            header = True
            print('Contig: ' + contig)
        if not contig in htable.keys():
            if (header):
                ofile.write(">C"+str(notplasmidcount)+'\n')
                notplasmidname[contig] = ("C"+str(notplasmidcount))
                notplasmidcount+=1
            else:
                ofile.write(line+'\n')
            header = False
        if contig in htable.keys():
            if (header):
                ofile.write(">" + htable.get(contig) + '\n')
            else:
                ofile.write(line+'\n')
            header = False

ofile.close()
print ("total contigs:",total)
print("writing table of paths: ",args.ofn_table)
ofile = open(args.ofn_table,"w+")
ofile.write("Plasmid\tIndex\tContig\tContig_len\tCumsum\n")
with open(args.ifn_recycler) as ifile:
    plasmidnum = 1
    for line in ifile:
        if (line.__contains__('RNODE')):
            plasmidnum = re.sub('_length.+', "", line)
            plasmidnum = int(plasmidnum[6:])
        if (line[0]=='(' or line[0] == "N"):
            line = line.rstrip("\n")
            line = line.replace("(", "").replace(",", "").replace("'", "").replace(")", "")
            line = line.replace("\"", "")
            contigs = line.split(" ")
            index = 1
            cumsum = 0
            for contig in contigs:
                print(contig)
                contig = contig[4:]
                contig = re.sub('^_[0-9]+_', "", contig)
                contig = re.sub('_ID_[0-9]+$', "", contig)
                print(htable.get(contig))
                print(plasmidnum)
                length = int(re.findall('^[0-9]+',contig[7:])[0])
                ofile.write("P"+str(plasmidnum)+'\t'+str(index)+'\t'+htable.get(contig)+'\t'+str(length)+'\t'+str(cumsum)+'\n')
                cumsum += length-k
                index+=1

# print(htable.keys())
# ofastg = open(args.ofn_fastg,"w+")
# with open(args.ifn_fastg) as fastg:
#     write = False
#     for line in fastg:
#         line = line.rstrip("\n")
#         if line[0] == '>':
#             apostrophe = []
#             line = line[1:]
#             line = line.replace(";", "")
#             line = line.split(':')
#             for lines in line:
#                 if lines.__contains__(','):
#                     line.remove(lines)
#                     lines = lines.split(',')
#                     for element in lines:
#                         line.append(element)
#             write = False
#             for contig in line:
#                 contig = contig[4:]
#                 if contig[-1]=="'":
#                     apostrophe.append("'")
#                 else:
#                     apostrophe.append("")
#                 contig = contig.replace("'","")
#                 contig = re.sub('^_[0-9]+_', "", contig)
#                 if contig in htable.keys():
#                     write = True
#         if write:
#             if line[0][0] == 'E':
#                 contignames = []
#                 contigcount = 0
#                 for contig in line:
#                     contig = contig[4:]
#                     contig = re.sub('^_[0-9]+_', "", contig)
#                     contig = contig.replace("'","")
#                     if contig in htable.keys():
#                         contignames.append(htable.get(contig))
#                     if contig in notplasmidname.keys():
#                         contignames.append(notplasmidname.get(contig))
#                     contigcount+=1
#                 ofastg.write('>')
#                 fastgcount = 1
#                 for contig in contignames:
#                     label = re.sub('_','-',contignames[fastgcount-1])
#                     apos=apostrophe[fastgcount-1]
#                     if fastgcount == contigcount:
#                         ofastg.write('EDGE_'+label+'_length_0_cov_0'+apos+";")
#                     elif fastgcount == 1:
#                         ofastg.write('EDGE_'+label+'_lenth_0_cov_0'+apos+":")
#                     else:
#                         ofastg.write('EDGE_'+label+'_length_0_cov_0'+apos+",")
#                     fastgcount+=1
#                 ofastg.write('\n')
#             else:
#                 ofastg.write(line+'\n')

