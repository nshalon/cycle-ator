import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Calculate read length distribtion and calculate coverage')
parser.add_argument('fasta', metavar='<tail_size>', type=str, help='percentile size of tail for insert distribution')
parser.add_argument('out_file', metavar='<paired_table>', type=str, help='paired table of reads')
args = parser.parse_args()

ofile = open(args.out_file,'w+')
ofile.write('cycle\tbp\tgc\n')
with open(args.fasta) as fasta:
    for line in fasta:
        if line[0] == '>':
            title = line.split(':')[0][5:]
        else:
            line = list(line)
            bp_count = 0
            gc_array = []
            gc_bp = []
            gc_count = 0
            length = len(line)
            bucket = int(length/250)
            if bucket < 10:
                bucket = 10
            for bp in line:
                bp_count += 1
                if bp == 'C' or bp == 'G':
                    gc_count += 1
                if bp_count % bucket == 0 or bp_count == length:
                    if bp_count == length and bp_count % bucket != 0:
                        gc_count = gc_count * (float(bucket)/(bp_count % bucket))
                        ofile.write(title + '\t' + str(bp_count) + '\t' + str(gc_count/float(bucket)) + '\n')
                        break
                    ofile.write(title + '\t' + str(bp_count) + '\t' + str(gc_count/float(bucket)) + '\n')
#                     gc_array.append(gc_count/bucket)
                    gc_count = 0
#                     gc_bp.append(bp_count)
                    
           #  graph = plt.plot(gc_bp,gc_array,'g')
#             plt.title(title)
#             plt.xlabel('bp')
#             plt.ylabel('G/C')
#             axes = plt.gca()
#             axes.set_ylim([0,1])
#             ofn_graph = args.out_dir + '/gc_' + title + '.pdf'
#             plt.savefig(ofn_graph)
#             plt.close()
                