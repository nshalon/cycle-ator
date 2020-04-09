import argparse
import re
#####################################################################################
# command line arguments
#####################################################################################

# define arguments
parser = argparse.ArgumentParser(description='remove long PCE')
parser.add_argument('ifn_recycler', metavar='<recycler ifn>', type=str, help='Input recycler file')
parser.add_argument('ofn_recycler', metavar='<recycler ofn>', type=str, help='Output recycler file')
args = parser.parse_args()

ofile = open(args.ofn_recycler,"w+")
with open(args.ifn_recycler) as ifile:
    max = 0
    print_to_ofile = True
    for line in ifile:
        if line[0] == "R":
            cyc_size = re.sub('RNODE_[0-9]+_length_','',line)
            cyc_size = int(re.sub('_.+','',cyc_size))
            if cyc_size>200000:
                print_to_ofile = False
                line_count = 0
        if print_to_ofile:
            ofile.write(line)
        if not print_to_ofile:
            line_count += 1
            if line_count == 4:
                print_to_ofile = True