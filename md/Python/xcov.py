import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math

parser = argparse.ArgumentParser(description='Calculate read length distribtion and calculate coverage')
parser.add_argument('tail_size', metavar='<tail_size>', type=float, help='percentile size of tail for insert distribution')
parser.add_argument('ifn_pairedtable', metavar='<paired_table>', type=str, help='paired table of reads')
parser.add_argument('ifn_tablepath', metavar='<table with path ifn>', type=str, help='Input table w/ plasmid contigs')
parser.add_argument('ofn_graph', metavar='<read length graph>', type=str, help='output directory (and filename) for graph')
parser.add_argument('ofn_paired', metavar = '<ofn_paired_dir>', type=str, help='Output classified paired table')
parser.add_argument('ofn_covstats', metavar = '<ofn_covstats>', type=str, help='Output coverage and read length stats')
parser.add_argument('k', metavar='<k>', type=int, help='kmer size used in assembly')

args = parser.parse_args()
max_perc = 100 - args.tail_size
min_perc = args.tail_size

print ("ifn paired table: ", args.ifn_pairedtable)
print ("ifn tablepath: ", args.ifn_tablepath)
print ("ofn read length graph: ", args.ofn_graph)
print ("ofn read length graph: ", args.ofn_paired)

print ("Top percentile: ", str(max_perc))
print ("Bottom percentile: ", str(min_perc))

print('Calculating insert sizes')

k = args.k

col_htable = {}
with open(args.ifn_pairedtable) as itable:
    line = itable.readline()
    header = line.split()
    index = 0
    for col in header:
        col_htable[col] = index
        index += 1
    
with open(args.ifn_pairedtable) as table:
    next(table)
    read_length_99 = read_length_1 = None
    insert_dist = np.array([],dtype=int)
    count = 0
    for line in table:
        count += 1 
        if count%1 == 0:
            line = line.split()
            if line[col_htable['contig1']] == line[col_htable['contig2']]:
                insert = abs(int(line[col_htable['back_coord1']])-int(line[col_htable['back_coord2']])) #switched to backcoords
                if (insert < 1600):    
                    insert_dist = np.append(insert_dist,insert)
    bins = []
    print(len(insert_dist))
    for n in range(np.amin(insert_dist),np.amax(insert_dist),3): bins.append(n)
    read_length_max = np.percentile(insert_dist, max_perc)
    read_length_min = np.percentile(insert_dist, min_perc)
    read_length_med = np.percentile(insert_dist, 50)
    hist = plt.hist(insert_dist, bins=bins)
    plt.title("Read Length Distribution ")
    plt.xlabel('Insert Size')
    plt.ylabel('Count')
    plt.axvline(x=read_length_max,color='r')
    plt.axvline(x=read_length_min,color='r')
    ofn_graph = args.ofn_graph + '/insert_dist.pdf'
    plt.savefig(ofn_graph)
    plt.close()
    table.close()

print(str(max_perc) + "% insert size:",int(read_length_max))
print(str(min_perc) + "% insert size:",int(read_length_min))

with open(args.ifn_tablepath) as path:
    for line in path: last=line
with open(args.ifn_tablepath) as path:
    contigtoplas = {}
    tables = []
    out_cyc_tables = []
    plaslength = 0
    plasmid = None 
    next(path)
    read_count = []
    read_count.append(0)
    for line in path:
        if line == last:
            last = True
        line = line.split()
        if (plasmid != line[0] and plasmid != None) or last==True:
            if(last == True and int(line[1])==1) and int(line[0][1:])==1:
                plaslength += int(line[3]) - k
                table = np.zeros(plaslength, dtype=int)  # create 2 2D tables, one for each strand of plasmid
                tables.append(np.copy(table))
                tables.append(np.copy(table))
                tables.append(np.copy(table))
                out_cyc_tables.append(np.copy(table))
                out_cyc_tables.append(np.copy(table))
                plaslength = 0
                print('1')
            if(last == True and int(line[1])==1):
                plaslength = int(line[3]) -k #REMOVE
                table = np.zeros(plaslength, dtype=int)  # create 2 2D tables, one for each strand of plasmid
                tables.append(np.copy(table))
                out_cyc_tables.append(np.copy(table))
                out_cyc_tables.append(np.copy(table))
                plaslength = 0
                print('2')
            if (last==True and int(line[0][1:])!=1) or (int(line[1])!=1):
                print('3')
                plaslength += int(line[3]) - k  # change the k value in both this script and the rename_path.py to accomodate for k = k
            table = np.zeros(plaslength, dtype=int)  # create 2 2D tables, one for each strand of plasmid
            tables.append(np.copy(table))
            out_cyc_tables.append(np.copy(table))
            out_cyc_tables.append(np.copy(table))
            print('Made table of size: ' + str(plaslength))
            plaslength = 0
            read_count.append(0)
        contigtoplas.setdefault(line[2], [])
        contigtoplas[line[2]].append(line[0])
        plasmid = line[0]
        plaslength += int(line[3]) - k
    path.close()

print(contigtoplas)

def cycle_to_plasmid(plasmid,contig,coord,strand):
    strand = int(strand)
    with open(args.ifn_tablepath) as plasmidpath:
        for line in plasmidpath:
            line = line.split()
            if line[2] == contig and line[0] == plasmid:
                start = int(line[4])
                orientation = line[5]
                length = int(line[3])
                break
        if (orientation == '-'):
            plascoord = start + (length - k - coord)
            if (strand == -1):
                t = 1
            else:
                t = -1
        else:
            plascoord = start + coord - k
            t = strand
        return [plascoord, t]

def get_start_end(plasmid,contig):
    with open(args.ifn_tablepath) as plasmidpath:
        for line in plasmidpath:
            line = line.split()
            if line[2] == contig and line[0] == plasmid:
                start = int(line[4])
                stop = start + int(line[3]) - k
                return [start,stop]

print("Now working on coverages...")
print("Length:",len(read_count))
with open(args.ifn_tablepath) as plasmidpath:
    out_cyc_count = 0
    cyc_count = 0
    max_dist = read_length_max
    min_dist = read_length_min
    next(plasmidpath)
    plasmid = None
    med_covs = []
    plasmid_num = -1
    for linepath in plasmidpath:
        linepath = linepath.split()
        if linepath[0] == '27':
            cyc_27 = True
        else:
            cyc_27 = False
        if(plasmid!=linepath[0] or plasmid==None):
            plasmid_num += 1
            plasmid = linepath[0]
            print("Now working on coverage for plasmid:",plasmid)
            with open(args.ifn_pairedtable) as table:
                next(table)
                for line in table:
                    line = line.split()
                    contig1 = line[1]
                    contig2 = line[14]
                    if plasmid in contigtoplas.get(contig1,"None") and plasmid in contigtoplas.get(contig2,"None"):
                        set1 = cycle_to_plasmid(plasmid,contig1,int(line[col_htable['back_coord1']]),int(line[col_htable['strand1']]))
                        set2 = cycle_to_plasmid(plasmid,contig1,int(line[col_htable['back_coord2']]),int(line[col_htable['strand2']]))
                        if set1[1] != set2[1]:
                            leftbackcoord = set1[0]
                            rightbackcoord = set2[0]
                            distance = (rightbackcoord - leftbackcoord) * set1[1]
                            min_coord = min(leftbackcoord , rightbackcoord)
                            max_coord = max(leftbackcoord , rightbackcoord)
                            if distance < 0:
                                distance = min_coord + len(tables[plasmid_num]) - max_coord
                                if distance < max_dist and distance > min_dist:
                                    read_count[plasmid_num] += 1 
                                    for bp in range(max_coord,len(tables[plasmid_num])):
                                        tables[plasmid_num][bp] += 1
                                        cyc_count += 1
                                    for bp in range(0,min_coord):
                                        cyc_count += 1
                                        tables[plasmid_num][bp] += 1
                            else:
                                if distance < max_dist and distance > min_dist:
                                    read_count[plasmid_num] += 1
                                    for bp in range(min_coord-1,max_coord-1):
                                        cyc_count += 1
                                        if cyc_27:
                                            print('here')
                                        if bp < len(tables[plasmid_num]):
                                            tables[plasmid_num][bp] += 1
                                        else:
                                            print('exceed plas count',plasmid,bp,len(tables[plasmid_num]))
                    elif (plasmid in contigtoplas.get(contig1,"None") and plasmid not in contigtoplas.get(contig2,"None")) or (plasmid not in contigtoplas.get(contig1,"None") and plasmid in contigtoplas.get(contig2,"None")):
                        if plasmid in contigtoplas.get(contig1,"None"):
                            set = cycle_to_plasmid(plasmid,contig1,int(line[col_htable['back_coord1']]),int(line[col_htable['strand1']]))
                            contig_in = contig1
                        else:
                            set = cycle_to_plasmid(plasmid,contig2,int(line[col_htable['back_coord2']]),int(line[col_htable['strand2']]))
                            contig_in = contig2   
                        start_and_stop_coords = get_start_end(plasmid,contig_in)
                        start = start_and_stop_coords[0]
                        stop = start_and_stop_coords[1]
                        if set[1] == 1:
                            minim = int(set[0])
                            if stop<(read_length_med+minim):
                                maxim = int(stop)
                            else:
                                maxim = int(read_length_med + minim)
                            leaving_reads_strand_num = 0
                        else:
                            maxim = int(set[0])
                            if start>(maxim-read_length_med):
                                minim = int(start)
                            else:
                                minim = int(maxim - read_length_med)
                            leaving_reads_strand_num = 1
                        for bp in range(minim-1,maxim):
                            out_cyc_count += 1
                            out_cyc_tables[2*plasmid_num+leaving_reads_strand_num][bp] += 1 
                            if cyc_27:
                                print('here')                                       
            med_covs.append(np.median(tables[plasmid_num]))
            print("Average coverage for plasmid:",int(np.median(tables[plasmid_num])),"SD:",int(np.std(tables[plasmid_num])))

with open(args.ofn_covstats,"w+") as stats:
    stats.write("PCE\tLength\tMed_cov\tNum_reads\tMed_read_dist\tMax_read_dist\tP_success\tCov_mean\tCov_SD\tBottleneck\tB_zscore\tB/M\tTwo_under\tThree_under\tFour_under\tFive_under\tSix_under\tSeven_under\n")
    for i in range(0,len(med_covs)):
        p_success = round((read_length_med/len(tables[i])),4)
        cov_mean = p_success*read_count[i]
        cov_sd = math.sqrt(med_covs[i])
        bottleneck = np.amin(tables[i])
        bn_z = round(((bottleneck-med_covs[i])/cov_sd),2)
        b_over_m = round(bottleneck/med_covs[i],2)
        two_under = med_covs[i] + (-2*cov_sd)
        three_under = med_covs[i] + (-3*cov_sd)
        four_under = med_covs[i] + (-4*cov_sd)
        five_under = med_covs[i] + (-5*cov_sd)
        six_under = med_covs[i] + (-6*cov_sd)
        seven_under = med_covs[i] + (-7*cov_sd)
        stats.write('P' + str(i+1) + '\t' + str(len(tables[i])) + '\t' + str(med_covs[i]) + '\t' + str(read_count[i]) + '\t' + str(read_length_med) + '\t' + str(round(read_length_max)) + '\t' + str(p_success) + '\t' + str(round(cov_mean,2)) + '\t' + str(round(cov_sd,2)) + '\t' + str(bottleneck) + '\t' + str(bn_z) + '\t' + str(b_over_m) + '\t' + str(two_under) + '\t' + str(three_under) + '\t' + str(four_under) + '\t' + str(five_under) + '\t' + str(six_under) + '\t' + str(seven_under) + '\n')

with open(args.ofn_paired,"w+") as results:
    plasmid_count = -1
    results.write('Plasmid\tBase_Pair\tCov\tOut_cyc_pos\tOut_cyc_neg\n')
    for i in range(0,len(med_covs)):
        plasmid_count += 1
        for bp in range(0,len(tables[plasmid_count])):
            cov = tables[plasmid_count][bp]
            out_cyc_pos = out_cyc_tables[2*plasmid_count][bp]
            out_cyc_neg = out_cyc_tables[2*plasmid_count+1][bp]
            if bp == (len(tables[plasmid_count])-1):
                out_cyc_pos = out_cyc_tables[2*plasmid_count][bp-1]
                out_cyc_neg = out_cyc_tables[2*plasmid_count+1][bp-1]
            if plasmid_count ==26:
                print(cov)
            results.write(str(plasmid_count+1)+'\t'+str(bp)+'\t'+str(cov)+'\t' + str(out_cyc_pos) + '\t' + str(out_cyc_neg) + '\n')
print("Find results @", args.ofn_paired)