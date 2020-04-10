all: plot_cycs

DIR_DONE=$(CYC_FIND_OUT)/.cyc_find_dir_done
$(DIR_DONE):
	$(_start)
	mkdir -p $(CYC_FIND_OUT)
	$(_end_touch)
dir: $(DIR_DONE)

ASSEMBLE_DONE?=$(CYC_FIND_OUT)/.assemble_done
$(ASSEMBLE_DONE): $(DIR_DONE)
	$(_start)
	$(_mega)/megahit -m 0.84 -o $(CYC_FIND_OUT)/megahit --bubble-level 0 --min-contig-len 27 --k-min 27 --k-max $(_k) --k-step 10 --merge-level 1000,0.95 -1 $(read1) -2 $(read2) -t 40
	$(_end_touch)
assemble: $(MF_DONE)

FASTA_TO_FASTG_DONE?=$(CYC_FIND_OUT)/.fasta_to_fastg
$(FASTA_TO_FASTG_DONE): $(ASSEMBLE_DONE)
	$(_start)
	$(_mega)/megahit_toolkit contig2fastg $(_k) $(CYC_FIND_OUT)/megahit/intermediate_contigs/k$(_k).contigs.fa > $(CYC_FIND_OUT)/k$(_k).fastg
	$(_end_touch)
fasta_to_fastg: $(FASTA_TO_FASTG_DONE)

PARSE_FASTG?=$(CYC_FIND_OUT)/.parse_fastg
$(PARSE_FASTG): $(FASTA_TO_FASTG_DONE)
	$(_start)
	python $(_py)/parse_fastg.py $(CYC_FIND_OUT)/k$(_k).fastg  $(CYC_FIND_OUT)/adjacency_list $(CYC_FIND_OUT)/contig_rename_map $(CYC_FIND_OUT)/renamed_final_contigs.fa $(CYC_FIND_OUT)/adjacency_matrix
	$(_end_touch)
parse_fastg: $(PARSE_FASTG)


GEN_K_SEQS?=$(CYC_FIND_OUT)/.kmer_fasta
$(GEN_K_SEQS): $(PARSE_FASTG)
	$(_start)
	python $(_py)/get_k_seqs.py $(CYC_FIND_OUT)/renamed_final_contigs.fa $(CYC_FIND_OUT)/adjacency_list $(_k) $(CYC_FIND_OUT)/kmer.fasta $(CYC_FIND_OUT)/kmer_list
	$(_end_touch)
gen_k_seqs: $(GEN_K_SEQS)

#MAP_K?=$(CYC_FIND_OUT)/.map_kmers
#$(MAP_K): #$(GEN_K_READS)
#	$(_start)
#	$(_bwa) index $(CYC_FIND_OUT)/mol.fasta
#	$(_bwa) mem $(CYC_FIND_OUT)/mol.fasta $(CYC_FIND_OUT)/kmer.fastq > $(CYC_FIND_OUT)/kmer-aln_ns.sam
#	perl $(_pl)/parse_bwa_sam.pl $(CYC_FIND_OUT)/kmer-aln_ns.sam $(CYC_FIND_OUT)/kmer_map_table_n $(CYC_FIND_OUT)/kmer_map_table_stats_n
#	#perl $(_pl)/filter_map.pl $(CYC_FIND_OUT)/kmer_map_table `expr $(_k) - 20` `expr $(_k) - 3` 3 $(CYC_FIND_OUT)/filter_kmer_map_table $(CYC_FIND_OUT)/filter_kmer_map_table_stats
#	$(_end_touch)
#map_k : $(MAP_K)

GET_WEIGHTS?=$(CYC_FIND_OUT)/.get_weight
$(GET_WEIGHTS): $(GEN_K_SEQS)
	$(_start)
	python $(_py)/precise_k_cov.py $(CYC_FIND_OUT)/contig_rename_map $(CYC_FIND_OUT)/kmer_list $(CYC_FIND_OUT)/k_reads.fastq $(CYC_FIND_OUT)/adjacency_list $(CYC_FIND_OUT)/kmer.fasta $(CYC_FIND_OUT)/edge_summary $(CYC_FIND_OUT)/edge_weight_matrix   
	$(_end_touch)
get_weights : $(GET_WEIGHTS) 

SHORT_CYCS?=$(CYC_FIND_OUT)/.short_cycs
$(SHORT_CYCS): $(GET_WEIGHTS)
	$(_start)
	mkdir -p $(CYC_FIND_OUT)/cycle_stats
	python $(_py)/get_shortest_cycs.py $(CYC_FIND_OUT)/adjacency_list $(CYC_FIND_OUT)/edge_summary $(CYC_FIND_OUT)/cycle_stats $(CYC_FIND_OUT)/summary_cycles
	cat $(CYC_FIND_OUT)/summary_cycles | cut -f1-6,12- | uniq > $(CYC_FIND_OUT)/summary_cycles_per_cyc
	$(_end_touch)
short_cycs: $(SHORT_CYCS)

PLOT_CYCS?=$(CYC_FIND_OUT)/.plot2
$(PLOT_CYCS): $(SHORT_CYCS)
	$(_start)
	mkdir -p $(CYC_FIND_OUT)/cycle_plots
	Rscript $(_r)/plot_putative_cycles.R $(CYC_FIND_OUT)/summary_cycles $(CYC_FIND_OUT)/cycle_stats $(CYC_FIND_OUT)/cycle_plots $(CYC_FIND_OUT)/adjacency_matrix $(CYC_FIND_OUT)/edge_summary 
	Rscript $(_r)/plot_cycle_coverages.R $(CYC_FIND_OUT)/summary_cycles $(CYC_FIND_OUT)
	$(_end_touch)
plot_cycs: $(PLOT_CYCS) 















#GEN_K_READS?=$(CYC_FIND_OUT)/.k_reads
#$(GEN_K_READS): $(PARSE_FASTG)
#       $(_start)
#       #       cat $(read1) $(read2) > $(CYC_FIND_OUT)/unpaired_reads.fastq
#       python $(_py)/generate_k_reads.py $(CYC_FIND_OUT)/unpaired_reads.fastq $(_k) $(CYC_FIND_OUT)/k_reads.fastq
#       #       $(_end_touch)
#gen_k_reads: $(GEN_K_READS)

#INDEX_FASTA?=$(CYC_FIND_OUT)/.index
#$(INDEX_FASTA): $(RENAME_FASTA)
#       $(_start)
#       #       $(_bwa) index $(CYC_FIND_OUT)/renamed_final_contigs.fa
#       $(_end_touch)
#       #index_fasta: $(INDEX_FASTA)
#
#MAP_DONE?=$(CYC_FIND_OUT)/.map_contigs
#$(MAP_DONE): $(INDEX_FASTA)
#       $(_start)
#       #       perl $(_pl)/trim_fastq.pl $(read1) 0 40 $(CYC_FIND_OUT)/trimR1.fastq
#       perl $(_pl)/trim_fastq.pl $(read2) 0 40 $(CYC_FIND_OUT)/trimR2.fastq
#       #       $(_bwa) mem $(CYC_FIND_OUT)/renamed_final_contigs.fa $(CYC_FIND_OUT)/trimR1.fastq > $(CYC_FIND_OUT)/r1-aln.sam
#       $(_bwa) mem $(CYC_FIND_OUT)/renamed_final_contigs.fa $(CYC_FIND_OUT)/trimR2.fastq > $(CYC_FIND_OUT)/r2-aln.sam
#       #       perl $(_pl)/parse_bwa_sam.pl $(CYC_FIND_OUT)/r1-aln.sam $(CYC_FIND_OUT)/R1table $(CYC_FIND_OUT)/R1table_stats
#       perl $(_pl)/parse_bwa_sam.pl $(CYC_FIND_OUT)/r2-aln.sam $(CYC_FIND_OUT)/R2table $(CYC_FIND_OUT)/R2table_stats
#       #       perl $(_pl)/filter_map.pl $(CYC_FIND_OUT)/R1table 30 40 2 $(CYC_FIND_OUT)/filter_R1table $(CYC_FIND_OUT)/filter_R1table_stats
#       perl $(_pl)/filter_map.pl $(CYC_FIND_OUT)/R2table 30 40 2 $(CYC_FIND_OUT)/filter_R2table $(CYC_FIND_OUT)/filter_R2table_stats
#       #       perl $(_pl)/pair_reads.pl $(CYC_FIND_OUT)/filter_R1table $(CYC_FIND_OUT)/filter_R2table $(CYC_FIND_OUT)/paired_table $(CYC_FIND_OUT)/paired_table_stats
#       $(_end_touch)
#       #map_done: $(MAP_DONE)
#
##CONTIG_COV?=$(CYC_FIND_OUT)/.contig_cov
#$(CONTIG_COV): $(MAP_DONE)
#       $(_start)
#       #       mkdir -p $(CYC_FIND_OUT)/contig_covs
#       python $(_py)/get_coverage.py $(CYC_FIND_OUT)/paired_table $(CYC_FIND_OUT)/contig_rename_map $(CYC_FIND_OUT)/contig_coverage $(CYC_FIND_OUT)/contig_covs
#       #       $(_end_touch)
#contig_coverage: $(CONTIG_COV)
