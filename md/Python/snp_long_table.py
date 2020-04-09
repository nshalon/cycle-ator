import argparse
import numpy as np

parser = argparse.ArgumentParser(description='remove reads that map to contigs that arent part of cycles')
parser.add_argument('med_cov_table', metavar='<med cov table for subjects>', type=str, help='subjects')
parser.add_argument('summary', metavar='<summary table>', type=str, help='for getting the cycles and the lengths')
parser.add_argument('base_in_dir', metavar='<base in dir>', type=str, help='directory to go into each snp table')
parser.add_argument('out_dir', metavar='<where to save all snp tables>', type=str, help='')
parser.add_argument('out_subjects', metavar='<file with all subject names>', type=str, help='')
args = parser.parse_args()

good_subjects = []
with open(args.med_cov_table) as subjects:
    all_subs = subjects.readline()
    all_subs = all_subs.split()
    for sub_count in range(0,len(all_subs)-1):
        sub_parts = all_subs[sub_count].split('_')
        time = sub_parts[1][-1]
        if time == '1':
            sub_partt2 = all_subs[sub_count+1].split('_')
            time2 = sub_partt2[1][-1]
            if time2 == '2':
                good_subjects.append(sub_parts[0])

with open(args.summary) as summ:
    next(summ)
    lengths = []
    for line in summ:
        line = line.split()
        lengths.append(line[2])

subfile = open(args.out_subjects,'w+')
bases = ['a','c','g','t']  
for sub in good_subjects:
    subfile.write(sub + '\n')
    print('Making table for subject',sub)
    path_to_out_file = args.out_dir + '/' + sub
    ofile = open(path_to_out_file,'w+')
    ofile.write('cycle\tcoord\tt1\tt2\tbase\n')
    cycle_count = 0
    for length in lengths:
        cycle_count += 1
        print(cycle_count)
        cyc_t1 = None
        cyc_t2 = None
        cycles_by_time = []
        cyc_t1 = np.zeros((int(length), 4), dtype=float)
        cyc_t2 = np.zeros((int(length), 4), dtype=float)
        cycles_by_time.append(cyc_t1)
        cycles_by_time.append(cyc_t2)
        for time in range(1,3):
            cyc_tables_ind = time-1
            subject_id = sub + '_t' + str(time)
            path_to_snp_table = args.base_in_dir + '/' + subject_id + '/' + 'snp_out/snp_full.tab'
            with open(path_to_snp_table) as snp:
                next(snp)
                for lines in snp:
                    line = lines.split()
                    header = line[0].split(':')
                    cycle = int(header[0][4:])
                    coord = int(line[1])-1
                    if cycle == cycle_count:
                        count_a = int(line[2])
                        count_c = int(line[3])
                        count_g = int(line[4])
                        count_t = int(line[5])
                        sum = float(count_a + count_c + count_g + count_t)
                        perc_a = round(count_a/sum,3)
                        perc_c = round(count_c/sum,3)
                        perc_g = round(count_g/sum,3)
                        perc_t = round(count_t/sum,3)
                        if perc_a != 0:
                            cycles_by_time[cyc_tables_ind][coord][0] = perc_a
#                             print(cycles_by_time[cyc_tables_ind][coord][0])
#                             print(perc_a)
                        if perc_c != 0:
                            cycles_by_time[cyc_tables_ind][coord][1] = perc_c
#                             print(cycles_by_time[cyc_tables_ind][coord][1])
#                             print(perc_c)
                        if perc_g != 0:
                            cycles_by_time[cyc_tables_ind][coord][2] = perc_g
#                             print(cycles_by_time[cyc_tables_ind][coord][2])
#                             print(perc_g)
                        if perc_t != 0:
                            cycles_by_time[cyc_tables_ind][coord][3] = perc_t
#                             print(cycles_by_time[cyc_tables_ind][coord][3])
#                             print(perc_t)
        for bp_num in range(0,int(length)):
            for bp in range(0,4):
                base = bases[bp]
                percents = [0,0]
                for cycle_ind in range(0,2):
                    percents[cycle_ind] = cycles_by_time[cycle_ind][bp_num][bp]
#                     print(cycles_by_time[cycle_ind][bp_num][bp])
                write_entry = False
                for percent in percents:
                    if percent != 0:
                        # print(percent)
                        write_entry = True
                if write_entry:   
                    ofile.write(str(cycle_count) + '\t' + str(bp_num) + '\t' + str(percents[0]) + '\t' + str(percents[1]) + '\t' + base + '\n')