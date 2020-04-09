import argparse

#############################################################################
# Purpose: input a fasta and parsed fastg graph file, output sequences of k overlaps between different contigs
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
#remember that it's a k-1 overlap, with k being the size used for assembly
k = args.k - 1


complement_bases = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
def get_reverse_complement(str):
    reverse_complement = ''
    for letter in str[::-1]:
        reverse_complement += complement_bases[letter]
    return reverse_complement

def get_comp_contig(contig_label):
    contig_dir = contig_label[-1]
    contig_strand = contig_label[-2]
    contig_name = contig_label[:-2]
    if contig_dir == '-':
        comp_dir = '+'
    elif contig_dir == '+':
        comp_dir = '-'
    if contig_strand == 'a':
        comp_strand = 'b'
    elif contig_strand == 'b':
        comp_strand = 'a'
    return contig_name + comp_strand + comp_dir
    
#put all the contigs with their names and seqs into a hashtable
contig_seq = {}
with open(args.ifn_fasta) as contigs:
    for contig_line in contigs:
        contig_line = contig_line.rstrip()
        if contig_line[0] == ">":
            contig = contig_line[1:]
        else:
            seq = contig_line
            contig_seq[contig] = seq


ofasta = open(args.ofn_fasta,'w+')
oklist = open(args.ofn_kmer_usage,'w+')
k_count = 0
node_to_k = {}
adj_line_count = 0
edge_to_k = {}
seq_to_k = {}
with open(args.contig_adjacency_list) as adj_list:
    for edges in adj_list:
        
        edges = edges.rstrip()
        edge = edges.split(':')
        parent = edge[0]
        childs = edge[1]
        parent_dir = parent[-1]
        #essentially if your contig can be the parent to a
        #variable edge
        
        if parent_dir == '+' and len(childs) > 2:
            contig_sequence = contig_seq.get(parent[:-2])
            if parent[-2] == 'b': #reverse complement
                #the whole kmer is the of length k+1 (k=original assembly k in eqn)
                #k_parent_seq is the first k bases of the kmer
                #the last base of the kmer belongs to the child seq
                k_parent_seq_rev_comp = contig_sequence[:(k+2)]
                k_parent_seq = get_reverse_complement(k_parent_seq_rev_comp)
            else:
                k_parent_seq = contig_sequence[-1*(k+2):] 
            comp_parent = get_comp_contig(parent)
            children = childs.split(',')
            for child in children: #to determine if this is a unique k or the same sequence as something previous
                
                forw_edge = (parent,child ) 
                comp_child = get_comp_contig(child)
                comp_edge = (comp_child,comp_parent)
                if edge_to_k.get(forw_edge) == None: 
                    
                    #if this edge hasn't been visited by reverse complement
                    child_contig_seq = contig_seq.get(child[:-2])
                    if child[-2] == 'b':
                        last_kmer_base = get_reverse_complement( child_contig_seq[ -1*k - 2 ] )
                    else:
                        last_kmer_base = child_contig_seq[k+1]
                    k_seq = k_parent_seq + last_kmer_base
                    comp_k_seq = get_reverse_complement(k_seq)
                    #check if the sequence has been written already   
                    if (seq_to_k.get(k_seq) != None):
                        k_label = seq_to_k[ k_seq ]
                    elif (seq_to_k.get(comp_k_seq) != None):
                        k_label = seq_to_k[ comp_k_seq ]
                    else:
                        k_count += 1
                        k_label = "k" + str(k_count)
                        seq_to_k[k_seq] = k_label
                        ofasta.write(">" + k_label + "\n" + k_seq + "\n") # write the fasta seq only for new kmers
                    edge_to_k[forw_edge] = k_label
                    edge_to_k[comp_edge] = k_label
                    
                    oklist.write(k_label + ":" + parent + ":" + child + "\n")
                    oklist.write(k_label + ":" + comp_child + ":" + comp_parent + "\n")

