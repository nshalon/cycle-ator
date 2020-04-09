import argparse

#define arguments
parser = argparse.ArgumentParser(description='remove reads that map to contigs that arent part of cycles')
parser.add_argument('in_ids', metavar='<ifn_paired table reads', type=str, help='Input paired table of reads')
# parser.add_argument('ids_in_one', metavar='<ifn_paired table reads', type=str, help='Input paired table of reads')
parser.add_argument('in_fq', metavar='<ifn_paired table reads', type=str, help='Input paired table of reads')
parser.add_argument('out_fq', metavar='<ifn_paired table reads w/ reads removed>', type=str, help='Output table')

args = parser.parse_args()

# id_table_one = []
# with open(args.ids_in_one) as read1:
#     for id in read1:
#         id = id.rstrip('\n')
#         id_table_one.append(id)

id_table = {}
with open(args.in_ids) as ids:
    count = 0
    for id in ids:
        id = id.rstrip('\n')
        # if id in id_table_one:
        id_table[id] = 1

ofile = open(args.out_fq,'w+')    
with open(args.in_fq) as in_fq:
    line_count = 0
    id_index = 0
    print_line = False
    count = 0
    for line in in_fq:
        if line_count%4 == 0:
            id = line.split()[0][1:]
            if id_table.get(id,0) == 1:
                print_line = True
                count += 1               
                print(count)
        if print_line:
            ofile.write(line)
        if line_count%4 == 3:
            print_line = False
        line_count += 1
        
            
            
        
        
        
