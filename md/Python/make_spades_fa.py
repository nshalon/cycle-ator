import argparse

parser = argparse.ArgumentParser(description='make spades like fasta')
parser.add_argument('in_mega', metavar='<in directory>', type=str, help='dir')
parser.add_argument('out_spades', metavar='<>', type=str, help='')
args = parser.parse_args()

ospades = open(args.out_spades,'w+')
with open(args.in_mega) as mega:
    node_count = 1
    for line in mega:
        if line[0] == '>':
            line = line.split()
            len = line[3][4:]
            header = '>NODE_' + str(node_count) + '_length_' + len + '_cov_100'
            ospades.write(header + '\n')
            node_count += 1
        else:
            ospades.write(line)            