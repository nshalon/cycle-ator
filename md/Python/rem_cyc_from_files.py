import argparse

parser = argparse.ArgumentParser(description='remove cycles from tablepath_ori and summary')
parser.add_argument('in_removedcycles', metavar='<in path>', type=str, help='dir')
parser.add_argument('in_path', metavar='<in path>', type=str, help='dir')
parser.add_argument('in_summ', metavar='<in summ>', type=str, help='dir')
parser.add_argument('out_path', metavar='<out>', type=str, help='')
parser.add_argument('out_summ', metavar='<out>', type=str, help='')
args = parser.parse_args()

removed_cycles = []
with open(args.in_removedcycles) as rem_cyc:
    for line in rem_cyc:
        line = line.rstrip()
        removed_cycles.append(line)

opath = open(args.out_path,'w+')
with open(args.in_path) as path:
    for line in path:
        plas = line.split()[0][1:]
        if plas not in removed_cycles:
            opath.write(line)
            
osumm = open(args.out_summ,'w+')
with open(args.in_summ) as summ:
    for line in summ:
        plas = line.split()[0][1:]
        if plas not in removed_cycles:
            osumm.write(line)

    