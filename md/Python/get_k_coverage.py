import argparse

#############################################################################
# Purpose: input a fasta and parsed fastg graph file
# Intput: size of k for assembly, a fasta of contigs, and an adjacency list parsed from the fastg
# Output: a fasta of sequences of kmers, with right name for kmer heading, to map against and a list
# of edges and the kmers that correspond to it.
# the output of fasta for kmers of overlapping segments 
#############################################################################

parser = argparse.ArgumentParser(description='rename reads to have FW and RV')
parser.add_argument('ifn_fasta', metavar='<contig fasta to be parsed for sequences>', type=str, help='input renamed fasta')
parser.add_argument('contig_adjacency_list', metavar='<adjacency list of contigs parsed from the fastg>', type=str, help='adjacency list of edges from assembly graph')
parser.add_argument('k', metavar='<size of k for assembly, k-1 overlaps between adjacent contigs>', type=int, help='input k integer used to assemble genome')
parser.add_argument('ofn_fasta', metavar='<fasta file with all the renamed k-1mers>', type=str, help='fasta with the k-1mers renamed to map against ')
parser.add_argument('ofn_kmer_usage', metavar='<output of k-1mers for each edge>', type=str, help='edges connected to which kmers they use')

args = parser.parse_args()

#rename inputs
contigs_fasta = args.ifn_fasta 
adj_list = args.contig_adjacency_list
#remember that it's a k-1 overlap, with k being the size used for assembly
k = args.k - 1


complement_bases = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
def get_reverse_complement(str):
    reverse_complement = ''
    for letter in str[::-1]:
        reverse_complement += complement_bases[letter]
    return reverse_complement
    
#put all the contigs with their names and seqs into a hashtable
contig_seq = {}
with open(contigs_fasta) as contigs:
    for contig_line in contigs:
        print(contig_line)
        contig_line = contig_line.rstrip()
        if contig_line[0] == ">":
            contig = contig_line[1:]
        else:
            seq = contig_line
            contig_seq[contig] = seq

print(contig_seq.keys())

ofasta = open(args.ofn_fasta,'w+')
oklist = open(args.ofn_kmer_usage,'w+')
k_count = 1
with open(adj_list) as adj_list:
    for edges in adj_list:
        edges = edges.rstrip()
        edge = edges.split(':')
        parent = edge[0]
        childs = edge[1]
        parent_dir = parent[-1]
        #essentially if your contig can be the parent to a
        #variable edge
        if parent_dir == '+' and len(childs) > 2:
            parent_strand = parent[-2]
            contig_sequence = contig_seq[parent[:-2]]
            k_seq = ''
            # return the reverse complement of the k-mer 
            # if the strand is b and direction is +
            # doesn't really matter if mapping, but just to keep
            # it consistent
            if parent_strand == "a":
                k_seq = contig_sequence[:k]
            else:
                k_seq = get_reverse_complement(contig_sequence[-1 * k:])
            k_label = 'k' + str(k_count)
            k_count += 1
            ofasta.write('>' + k_label + '\n' + k_seq + '\n')
            for child in childs.split(','):
                oklist.write(k_label + ":" + parent + ":" + child + "\n")
            
        
    









