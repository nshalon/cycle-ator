import argparse

parser = argparse.ArgumentParser(description='remove reads that map to contigs that arent part of cycles')
parser.add_argument('in_snp', metavar='<in directory>', type=str, help='dir')
parser.add_argument('cycles', metavar='<amount of cycles>', type=int, help='amount of cycles to do this for')
parser.add_argument('out_snp', metavar='<>', type=str, help='')
args = parser.parse_args()

print("snp table:",args.in_snp)

cycle_count = args.cycles + 1
not_ref = [True,True,True,True]
osnp = open(args.out_snp,'w+')
osnp.write('Cycle\tcoord\tA\tC\tG\tT\n')
for cycle_count in range(1,cycle_count):
    with open(args.in_snp) as snp:
        next(snp)
        for lines in snp:
            line = lines.split()
            header = line[0].split(':')
            cycle = int(header[0][4:])
            coord = line[1]
            if cycle == cycle_count:
                ref = line[6]
                count_a = int(line[2])
                count_c = int(line[3])
                count_g = int(line[4])
                count_t = int(line[5])
                sum = count_a + count_c + count_g + count_t
                max_morph = max(count_a,count_c,count_g,count_t)
                if sum<=10:
                    ratio = 0.7
                elif sum>100:
                    ratio = 0.95
                else:
                    ratio = 0.9
                if max_morph/sum < ratio and sum>5:
                    osnp.write(str(cycle)+'\t'+coord+'\t' + str(count_a) + '\t' + str(count_c) + '\t' + str(count_g) + '\t' + str(count_t) + '\n')               
                else:
                    print(lines)
            