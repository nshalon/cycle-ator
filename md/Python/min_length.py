import argparse

#############################################################################
# Purpose: find shortest cycle back to (-) node of each (+) node. Start point is (+) node for each (a) strand
# Input: all nodes parsed and renamed from fastg, edges, and weights
# Output: table of all shortest cycles, a separate adjacency matrix for each, and a weight list for r graphs
#############################################################################

parser = argparse.ArgumentParser(description='rename reads to have FW and RV')
parser.add_argument('ifn_edges', metavar='<list of all renamed edges from fastg>', type=str, help='input fastq files')


args = parser.parse_args()

min_length = 10000000000
linec = 0
with open(args.ifn_edges) as edges:
    for line in edges:
        if line[0] == ">":
            linec += 1
            if linec % 100000 == 0:
                print(linec)
            line = line.split()
            length_str = line[3]
            length = length_str[4:]
            if int(length) < min_length:
                min_length = int(length)
print("MIN LENGTH:",min_length)
                