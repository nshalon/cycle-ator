#!/usr/bin/env python

import argparse
import numpy as np

#define arguments
parser = argparse.ArgumentParser(description='Map reads to plasmids and classify')
parser.add_argument('ifn_pairedtable', metavar='<pairedtable ifn>', type=str, help='Input paired table')
parser.add_argument('ifn_tablepath', metavar='<table with path ifn>', type=str, help='Input table w/ plasmid contigs')
parser.add_argument('distancecutoff', metavar='<distance_cut_off>', type=int, help='The maximum distance bw paired reads')
parser.add_argument('ofn_table', metavar='<ofn_table>', type=str, help='Output tables')
parser.add_argument('ofn_paired', metavar = '<ofn_paired_dir', type=str, help='Output classified paired table')
parser.add_argument('k', metavar='<k>', type=int, help='kmer size used in assembly')
parser.add_argument('stats', metavar='<stats>', type=str, help='far reads debug file')
#parse arguments
args = parser.parse_args()
print ("ifn paired table: ", args.ifn_pairedtable)
print ("ifn tablepath: ", args.ifn_tablepath)
print ("distance_cut_off: ", args.distancecutoff)
print ("Output table: ", args.ofn_table)
print("Output paired table directory: ", args.ofn_paired)

print("Reading the table w/ pasmid contigs file: ", args.ifn_tablepath)

k = args.k

col_htable = {}
with open(args.ifn_pairedtable) as itable:
    line = itable.readline()
    header = line.split()
    index = 0
    for col in header:
        col_htable[col] = index
        index += 1
        
print(col_htable)

print(col_htable['contig1'],col_htable['contig2'])
print(col_htable['coord1'],col_htable['coord2'])
########################
#Create htable w/ keys as contigs and values as plasmids and creates list of tables the size of each plasmid
#################

with open(args.ifn_tablepath) as plasmidpath:
    num_lines = sum(1 for line in open(args.ifn_tablepath))-1
    contigtoplas = {}
    plastotablenum = {}
    tablenumtoplas = {}
    tablenum = 0
    tables = []
    plasmids = []
    next(plasmidpath)
    plaslength = 0
    plasmid = None
    linenum = 0
    for line in plasmidpath:
        line = line.split()
        linenum+=1
        print(linenum,num_lines)
        if (plasmid != line[0] and plasmid != None) or (linenum == num_lines):
            if(num_lines == 1):
                plasmid = line[0]
            if(linenum==num_lines and int(line[1]) != 1):
                plaslength += int(line[3])-k #change the k value in both this script and the rename_path.py to accomodate for k = k
            plasmids.append(plasmid)
            plastotablenum[plasmid] = tablenum
            print('Plasmid:'+plasmid)
            tablenumtoplas[tablenum] = plasmid
            tablenumtoplas[tablenum+1] = plasmid
            print(plaslength)
            for i in range(0,2):
                table = np.zeros((4, plaslength), dtype=int)  # create 2 2D tables, one for each strand of plasmid
                tables.append(table)
                print('made table:' + str(plaslength))
            plaslength = 0
            tablenum += 2
            if(linenum==num_lines and int(line[1]) == 1):
                plasmid = line[0]
                plaslength += int(line[3])-k
                plasmids.append(plasmid)
                plastotablenum[plasmid] = tablenum
                print('Plasmid:'+plasmid)
                tablenumtoplas[tablenum] = plasmid
                tablenumtoplas[tablenum+1] = plasmid
                print(plaslength)
                for i in range(0,2):
                    table = np.zeros((4, plaslength), dtype=int)  # create 2 2D tables, one for each strand of plasmid
                    tables.append(table)

        contigtoplas.setdefault(line[2], [])
        contigtoplas[line[2]].append(line[0])
        plasmid = line[0]
        plaslength += int(line[3])-k

print(contigtoplas)


def cycle_to_plasmid(plasmid,contig,coord,strand,one):
    strand = int(strand)
    with open(args.ifn_tablepath) as plasmidpath:
        for line in plasmidpath:
            line = line.split()
            if line[2] == contig and line[0] == plasmid:
                start = int(line[4])
                orientation = line[5]
                length = int(line[3])
                break
        if(orientation == '-'):
            plascoord = start + (length - coord)
            if coord<k and one:
                return 'NULL'
            if(strand == -1):
                t = 1
            else:
                t = -1
        else: 
            if length - coord < k and one:
                return 'NULL'    
            plascoord = start + coord
            t = strand
    return [plascoord,t]
# cycle_to_plasmid('R1')


# def plasmid_to_cycle(plasmid,coord,t):
#     with open(args.ifn_tablepath) as plasmidpath:
#         length=0
#         for line in plasmidpath:
#             line = line.split()
#             if line[0] == plasmid:
#                 length += int(line[3])
#                 if length > coord:
#                     cumsum = int(line[4])
#                     orientation = line[5]
#                     contig = line[2]
#                     contiglen = int(line[3])
#                     break
#         j = coord - cumsum + 30
#         if(orientation == '-'):
#             if t == -1:
#                 strand = 1
#             else:
#                 strand = -1
#             j = contiglen-j
#         else:
#             strand = t
#     return [contig,j,strand]


stats = open(args.stats,'w+')
stats.write('distance,coord1,coord2,contig1,contig2,plasstrand1,plasstrand2\n')
print("\nNow reading the paired_table and editing matrices...")
def maketables(index1,index2): #put in function so that both columns of table can easily be read and classified
    with open(args.ifn_tablepath) as plasmidpath:
        next(plasmidpath)
        maxdist = args.distancecutoff
        plasmid = None
        for linepath in plasmidpath:
            linepath = linepath.split()
            contig = linepath[2]
            if(plasmid!=linepath[0] or plasmid==None):
                plasmid = linepath[0]
                print("Now working on plasmid:",plasmid)
                ofile = open((args.ofn_paired+'/'+plasmid+'_pairedtable'),"w+")
                ofile.write("Key\tContig\tCoord\tStrand\tPlasmid\tPlasmid_coord\tPlasmid_strand\tPair_contig\tPair_coord\tPair_strand\tPair_plasmid\tPair_plasmid_coord\tPair_plasmid_strand\tClass\n")
            count = 0
            with open(args.ifn_pairedtable) as table:
                next(table)
                for line in table:
                    line = line.split()
                    contig1 = line[index1]
                    if contig1 == contig:
                        contig2 = line[index2]
                        coord1 = int(line[index1+1])
                        coord2 = int(line[index2+1])
                        strand1 = int(line[index1 + 3])
                        strand2 = int(line[index2 + 3])
                        set1 = cycle_to_plasmid(plasmid,contig1,coord1,strand1,True)
                        if set1 == 'NULL':
                            continue                          
                        contig2list = contigtoplas.get(contig2)  # to make sure contig is in a plasmid
                        if(contig2list == None):
                            set2 = [-1,-1,-1]
                            contig2list = ['None']
                            plasmid2='not_plas'
                        else:
                            if(plasmid in contig2list): #if it's a repeat seq that appears in first plasmid
                                plasmid2 = plasmid
                            else:                       #if the contig doesn't appear in first plasmid then the plasmid is just the first plasmid in the list (repeat and non-repeats)
                                plasmid2 = contig2list[0]
                            set2 = cycle_to_plasmid(plasmid2,contig2,coord2,strand2,False)                        
                        plascoord1 = set1[0]-1
                        plascoord2 = set2[0]-1
                        plasstrand1 = set1[1]
                        plasstrand2 = set2[1]
                        oppositestrand = plasstrand1*plasstrand2<0
                        plasmidnum = 2 * (int(plasmid[1:]) - 1)
                        strandnum = 0
                        if (plasstrand1) == -1:
                            strandnum = 1                    
                        tablenum = plastotablenum.get(plasmid) + strandnum
                        if plasmid in contig2list and oppositestrand:
                            # if plascoord1 > len(tables[tablenum][1]):
#                                 print('Here',contig1,coord1,tablenum,len(tables[tablenum][1]))
                            distance = (plascoord2 - plascoord1) * int(plasstrand1)
                            if(distance<0): #if the reads cross the beginning or end of cycle
                                if plasstrand1 == -1:
                                    plascoord1sub = plascoord1 + len(tables[tablenum][0])
                                else:
                                    plascoord1sub = plascoord1 - len(tables[tablenum][0])
                                distance = (plascoord2 - plascoord1sub) * int(plasstrand1)
                                # print(str(distance) + ' negative')
                           #  if plascoord1<0:
#                                 print(plascoord1,coord1)
                            if distance <= maxdist and distance > 0:                 #good reads
                                # print(plascoord1,len(tables[tablenum][0]))
                                tables[tablenum][0][plascoord1] += 1
                                label = 'Good_Reads'
                            elif distance > 0:                                   #for reads far from each other
                                tables[tablenum][1][plascoord1] += 1
                                label = 'Distant_Reads'
                                line_stats = str(distance)+'\t'+str(coord1)+'\t'+str(coord2)+'\t'+contig1+'\t'+contig2+'\t'+str(plasstrand1)+'\t'+str(plasstrand2)+'\n'
                                stats.write(line_stats)                              
                            else:
                                label = 'Weird read'
                        elif contig2list[0] != 'None':
                            tables[tablenum][2][plascoord1] += 1
                            label = 'Not_on_plasmid'
                        else:                                       #reads on diff plasmids or same strand
                            tables[tablenum][3][plascoord1] += 1
                            label = 'On_lin_seq'
                        ofile.write(line[0]+'\t'+contig1+'\t'+str(coord1)+'\t'+str(strand1)+'\t'+plasmid+'\t'+str(plascoord1)+'\t'+str(plasstrand1)+'\t'+contig2+'\t'+str(coord2)+'\t'+str(strand2)+'\t'+plasmid2+'\t'+str(plascoord2)+'\t'+str(plasstrand2)+'\t'+label+'\n')

maketables(col_htable['contig1'],col_htable['contig2']) #These two are so that both columns of paired table are read and classified
maketables(col_htable['contig2'],col_htable['contig1'])

print("Writing output table in ", args.ofn_table)
with open(args.ofn_table, "w+") as out:
    out.write('Plasmid\tBase_Pair\tStrand\tGood_Reads\tDistant_Reads\tNot_on_plasmid\tOn_lin_seq')
    tablecount = 0
    for table in tables:
        plasmidnum = (tablenumtoplas.get(tablecount))[1:]
        strand = '+'
        if(tablecount%2 != 0):
            strand = '-'
        for i in range(0,len(table[0])):
            bp = str(i + 1)
            one = str(table[0][i])
            two = str(table[1][i])
            three = str(table[2][i])
            four = str(table[3][i])
            out.write('\n'+plasmidnum+'\t'+bp+'\t'+strand+'\t'+one+'\t'+two+'\t'+three+'\t'+four)
        print("Plasmid", plasmidnum, strand, "finished")
        tablecount+=1
print("Find the results @", args.ofn_table)
