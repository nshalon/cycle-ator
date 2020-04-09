import argparse

#############################################################################
# Purpose: given a list of reads, output all reads broken up into k-mers in fastq format
# Intput: size of k for assembly, fastq of all reads. Reads will be treated as unpaired
# Output: fastq format, each read the size of k which was from the original read
#############################################################################

parser = argparse.ArgumentParser(description='rename reads to have FW and RV')
parser.add_argument('ifn_fastq', metavar='<fastq of all one sided reads>', type=str, help='input fastq files')
parser.add_argument('k', metavar='<size of k for assembly, k-1 overlaps between adjacent contigs>', type=int, help='input k integer used to assemble genome')
parser.add_argument('ofn_k_fastq', metavar='<output fastq of sliding k window >', type=str, help='fasta with the k-1mers renamed to map against ')

args = parser.parse_args()

k = args.k - 1
ifastq = args.ifn_fastq

ofastq = open(args.ofn_k_fastq,'w+')
line_count = 1
header = ''
quality = ''
with open(ifastq) as ifastq:
    for line in ifastq:
        line = line.rstrip()
        if line_count == 1:
            header = line
        if line_count == 2:
            sequence = line
        if line_count == 4:
            quality = line
            k_count = 1
            if (len(sequence) < k) or ('N' in sequence):
                line_count = 1
                continue
                
            for start_position in range(0,( len(sequence) - (k + 2) ) ):
                sequence_subset = sequence[start_position:start_position+k+3]
                quality_subset = quality[start_position:start_position+k+3]
                subheader = header + "_" + str(k_count)
                k_count += 1
                ofastq.write(subheader + '\n' +
                                sequence_subset + '\n' +
                                '+' + '\n' + 
                                quality_subset + '\n')
            line_count = 1
            continue
        line_count += 1
    
    