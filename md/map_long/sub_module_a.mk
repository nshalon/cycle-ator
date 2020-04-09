all: graph_cov2

INDEX?=$(FIG_OUT_DIR2)/'.index'
$(INDEX):
	$(_bwa) index $(FIG_OUT_DIR)/renamefasta
index: $(INDEX)

MAP_OLDDONE?=$(FIG_OUT_DIR2)/'map_done2'
$(MAP_OLDDONE):
	$(_start)
	mkdir -p $(FIG_OUT_DIR2)
	$(_bwa) index $(FIG_OUT_DIR)/renamefasta
	#gunzip -c -f $(2_read1) > $(FIG_OUT_DIR2)/R1.fastq
	#gunzip -c -f $(2_read2) > $(FIG_OUT_DIR2)/R2.fastq
	cp $(2_read1) $(FIG_OUT_DIR2)/R1.fastq
	cp $(2_read2) $(FIG_OUT_DIR2)/R2.fastq
	perl $(_pl)/trim_fastq.pl $(FIG_OUT_DIR2)/R1.fastq 0 40 $(FIG_OUT_DIR2)/trimR1.fastq
	perl $(_pl)/trim_fastq.pl $(FIG_OUT_DIR2)/R2.fastq 0 40 $(FIG_OUT_DIR2)/trimR2.fastq
	$(_bwa) mem -t 40 $(FIG_OUT_DIR)/renamefasta $(FIG_OUT_DIR2)/trimR1.fastq > $(FIG_OUT_DIR2)/r1-aln.sam
	$(_bwa) mem -t 40 $(FIG_OUT_DIR)/renamefasta $(FIG_OUT_DIR2)/trimR2.fastq > $(FIG_OUT_DIR2)/r2-aln.sam
	perl $(_pl)/parse_bwa_sam.pl $(FIG_OUT_DIR2)/r1-aln.sam $(FIG_OUT_DIR2)/R1table $(FIG_OUT_DIR2)/R1table_stats
	perl $(_pl)/parse_bwa_sam.pl $(FIG_OUT_DIR2)/r2-aln.sam $(FIG_OUT_DIR2)/R2table $(FIG_OUT_DIR2)/R2table_stats
	perl $(_pl)/filter_map.pl $(FIG_OUT_DIR2)/R1table 28 36 2 $(FIG_OUT_DIR2)/filter_R1table $(FIG_OUT_DIR2)/filter_R1table_stats
	perl $(_pl)/filter_map.pl $(FIG_OUT_DIR2)/R2table 28 36 2 $(FIG_OUT_DIR2)/filter_R2table $(FIG_OUT_DIR2)/filter_R2table_stats
	perl $(_pl)/pair_reads.pl $(FIG_OUT_DIR2)/filter_R1table $(FIG_OUT_DIR2)/filter_R2table $(FIG_OUT_DIR2)/paired_table_pre $(FIG_OUT_DIR2)/paired_table_stats
	$(_end_touch)
mapold: $(MAP_OLDDONE)

REM_NCYC?=$(FIG_OUT_DIR)/'.rem_done'
$(REM_NCYC): $(MAP_OLDDONE)
	$(_start)
	python $(_py)/rem_noncyc_reads.py $(FIG_OUT_DIR2)/paired_table_pre $(FIG_OUT_DIR2)/paired_table $(FIG_OUT_DIR2)/paired_table_remove_stats
	$(_end_touch)
remove_reads: $(REM_NCYC)


CLASS_DONE2?=$(FIG_OUT_DIR2)/'class_done2'
$(CLASS_DONE2): $(REM_NCYC)
	$(_start)
	mkdir -p $(FIG_OUT_DIR2)/pairedtables
	python $(_py)/map_contigs.py $(FIG_OUT_DIR2)/paired_table $(FIG_OUT_DIR)/tablepath_ori 2000 $(FIG_OUT_DIR2)/classified_reads $(FIG_OUT_DIR2)/pairedtables $(_k) $(FIG_OUT_DIR2)/read_stats
	$(_end_touch)
class2: $(CLASS_DONE2)

BUCK_DONE2?=$(FIG_OUT_DIR2)/'bucket_done2'
$(BUCK_DONE2): $(CLASS_DONE2)
	$(_start)
	python $(_py)/table_bucket.py $(FIG_OUT_DIR2)/classified_reads 100 $(FIG_OUT_DIR2)/bucket_classified_reads_100
	python $(_py)/table_bucket.py $(FIG_OUT_DIR2)/classified_reads 10 $(FIG_OUT_DIR2)/bucket_classified_reads_10
	$(_end_touch)
bucket2: $(BUCK_DONE2)

PLOT_DONE2?=$(FIG_OUT_DIR2)/'plot_done2'
$(PLOT_DONE2): $(BUCK_DONE2)
	$(_start)
	mkdir -p $(GRAPH_OUT_DIR2)
	mkdir -p $(GRAPH_OUT_DIR2)/bucket100
	mkdir -p $(GRAPH_OUT_DIR2)/bucket10
	Rscript $(_r)/plot_table.r $(FIG_OUT_DIR2)/bucket_classified_reads_100 $(FIG_OUT_DIR)/tablepath_ori 100 $(GRAPH_OUT_DIR2)/bucket100
	Rscript $(_r)/plot_table.r $(FIG_OUT_DIR2)/bucket_classified_reads_10 $(FIG_OUT_DIR)/tablepath_ori 10 $(GRAPH_OUT_DIR2)/bucket10
	$(_end_touch)
graph2: $(PLOT_DONE2)

COV_DONE2?=$(GRAPH_OUT_DIR2)/'cov_done2'
$(COV_DONE2): $(REM_NCYC)
	$(_start)
	mkdir -p $(GRAPH_OUT_DIR2)/cov
	python $(_py)/xcov.py 1 $(FIG_OUT_DIR2)/paired_table $(FIG_OUT_DIR)/tablepath_ori $(GRAPH_OUT_DIR2)/cov $(GRAPH_OUT_DIR2)/cov_table $(GRAPH_OUT_DIR2)/stats $(_k)
	$(_end_touch)
coverage2: $(COV_DONE2)

GRAPH_COV2?=$(GRAPH_OUT_DIR2)/'graph_done2'
$(GRAPH_COV2): $(COV_DONE2)
	$(_start)
	Rscript $(_r)/graph_cov.R $(GRAPH_OUT_DIR2)/cov_table 1 $(FIG_OUT_DIR)/tablepath_ori $(GRAPH_OUT_DIR2)/cov $(GRAPH_OUT_DIR2)/stats $(_k)
	$(_end_touch)
graph_cov2: $(GRAPH_COV2)

SNP_FILES?=$(FIG_OUT_DIR)/'.snp'
$(SNP_FILES): #$(GRAPH_COV2)
	$(_start)
	mkdir -p $(FIG_OUT_DIR)/snp_in
	python $(_py)/make_cyc_table.py $(OUT_DIR)/unfil_recycle/k$(_k).cycs.fasta $(FIG_OUT_DIR)/snp_in/cycle_table $(FIG_OUT_DIR)/snp_in/cycle.fa
	$(_end_touch)
snp_files: $(SNP_FILES)

MAP_SNP?=$(FIG_OUT_DIR2)/'.snp_map'
$(MAP_SNP): #$(SNP_FILES)
	$(_start)
	
	mkdir -p $(FIG_OUT_DIR2)/snp_in
	#$(_bwa) index $(FIG_OUT_DIR)/snp_in/cycle.fa
	perl $(_pl)/trim_fastq.pl $(FIG_OUT_DIR2)/R1.fastq 0 110 $(FIG_OUT_DIR2)/snp_in_trimR1
	#perl $(_pl)/trim_fastq.pl $(FIG_OUT_DIR2)/R2.fastq 0 40 $(FIG_OUT_DIR2)/snp_in_trimR2
	$(_bwa) mem -t 40 $(FIG_OUT_DIR)/snp_in/cycle.fa $(FIG_OUT_DIR2)/snp_in_trimR1 > $(FIG_OUT_DIR2)/snp_r1-aln.sam
	#$(_bwa) mem -t 40 $(FIG_OUT_DIR)/snp_in/cycle.fa $(FIG_OUT_DIR2)/snp_in_trimR2 > $(FIG_OUT_DIR2)/snp_r2-aln.sam
	perl $(_pl)/parse_bwa_sam.pl $(FIG_OUT_DIR2)/snp_r1-aln.sam $(FIG_OUT_DIR2)/snp_in/R1.fastq $(FIG_OUT_DIR2)/snp_in/snp_R1table_stats
	#perl $(_pl)/parse_bwa_sam.pl $(FIG_OUT_DIR2)/snp_r2-aln.sam $(FIG_OUT_DIR2)/snp_in/R2.fastq $(FIG_OUT_DIR2)/snp_in/snp_R2table_stats
	$(_end_touch)
map_snp: $(MAP_SNP)

RUN_SNP?=$(FIG_OUT_DIR2)/'.runsnp'
$(RUN_SNP): #$(MAP_SNP)
	$(_start)
	mkdir -p $(FIG_OUT_DIR2)/snp_out
	mkdir -p $(FIG_OUT_DIR2)/snp_out/out_full
	mkdir -p $(FIG_OUT_DIR2)/snp_out/out_clipped
	$(varisum_dir)/varisum -idir $(FIG_OUT_DIR2)/snp_in -contigs $(FIG_OUT_DIR)/snp_in/cycle_table -contig_field contig -contigs_fa $(FIG_OUT_DIR)/snp_in/cycle.fa -min_score 30 -max_edit 2 -min_length 40 -odir_full $(FIG_OUT_DIR2)/snp_out/out_full -odir_clipped $(FIG_OUT_DIR2)/snp_out/out_clipped -ofn_snp_full $(FIG_OUT_DIR2)/snp_out/snp_full.tab -ofn_snp_clipped $(FIG_OUT_DIR2)/snp_clipped.tab
	$(_end_touch)
run_snp: $(RUN_SNP)

PARSE_SNP?=$(FIG_OUT_DIR2)/'.snp'
$(PARSE_SNP): #$(RUN_SNP)
	$(_start)
	python $(_py)/parse_snp.py $(FIG_OUT_DIR2)/snp_out/snp_full.tab 195 $(FIG_OUT_DIR2)/parsed_snp.tab
	$(_end_touch)
parse_snp: $(PARSE_SNP)

GRAPH_SNP?=$(FIG_OUT_DIR2)/'.graph_pce_snp'
$(GRAPH_SNP):
	$(_start)
	mkdir -p $(GRAPH_OUT_DIR)/pce_snp
	mkdir -p $(GRAPH_OUT_DIR)/good_cycles
	mkdir -p $(GRAPH_OUT_DIR)/good_cycles/pce_snp
	Rscript $(_r)/plot_snp.r $(FIG_OUT_DIR)/tablepath_ori $(GRAPH_OUT_DIR)/pce_snp $(GRAPH_OUT_DIR)/med_cov_table $(FIG_OUT_DIR) $(GRAPH_OUT_DIR)/summary $(GRAPH_OUT_DIR)
	Rscript $(_r)/plot_snp.r $(FIG_OUT_DIR)/tablepath_ori_rem $(GRAPH_OUT_DIR)/good_cycles/pce_snp $(GRAPH_OUT_DIR)/med_cov_table $(FIG_OUT_DIR) $(GRAPH_OUT_DIR)/summary_removed $(GRAPH_OUT_DIR)
	$(_end_touch)
graph_snp: $(GRAPH_SNP)

MAKE_SNP_LONG?=$(FIG_OUT_DIR)/'.longtabldone'
$(MAKE_SNP_LONG):
	$(_start)
	mkdir -p $(FIG_OUT_DIR)/snp_long_tables
	python $(_py)/snp_long_table.py $(GRAPH_OUT_DIR)/med_cov_table $(GRAPH_OUT_DIR)/summary $(FIG_OUT_DIR) $(FIG_OUT_DIR)/snp_long_tables $(FIG_OUT_DIR)/snp_long_tables/subjects
	$(_end_touch)
make_snp_long: $(MAKE_SNP_LONG)

GRAPH_SNP_LONG?=$(FIG_OUT_DIR)/'.snps_long_graph'
$(GRAPH_SNP_LONG): #$(MAKE_SNP_LONG)
	$(_start)
	mkdir -p $(GRAPH_OUT_DIR)/snp_long
	Rscript $(_r)/graph_snp_long.r $(GRAPH_OUT_DIR)/summary $(FIG_OUT_DIR)/snp_long_tables/subjects $(FIG_OUT_DIR)/snp_long_tables $(GRAPH_OUT_DIR)/snp_long
	$(_end_touch)
graph_snp_long: $(GRAPH_SNP_LONG)

REM_CYC?=$(GRAPH_OUT_DIR)/'.rem_cyc'
$(REM_CYC):
	$(_start)
	python $(_py)/rem_cyc.py $(GRAPH_OUT_DIR) $(GRAPH_OUT_DIR)/removed_cycles
	python $(_py)/rem_cyc_from_files.py $(GRAPH_OUT_DIR)/removed_cycles $(FIG_OUT_DIR)/tablepath_ori $(GRAPH_OUT_DIR)/summary $(FIG_OUT_DIR)/tablepath_ori_rem $(GRAPH_OUT_DIR)/summary_removed
	$(_end_touch)
rem_cyc: $(REM_CYC)

REM_CYCLES_GRAPH?=$(GRAPH_OUT_DIR)/'.graphed_rem_cyc'
$(REM_CYCLES_GRAPH):
	$(_start)
	mkdir -p $(GRAPH_OUT_DIR)/good_cycles
	mkdir -p $(GRAPH_OUT_DIR)/good_cycles/PCE_all_subs
	Rscript $(_r)/graph_long.r $(num_sub) $(FIG_OUT_DIR)/tablepath_ori_rem $(GRAPH_OUT_DIR) $(GRAPH_OUT_DIR)/good_cycles/PCE_all_subs $(_k) $(GRAPH_OUT_DIR)/med_cov_table
	Rscript $(_r)/all_pce_and_subs.r $(FIG_OUT_DIR)/tablepath_ori_rem $(GRAPH_OUT_DIR) $(GRAPH_OUT_DIR)/med_cov_table 
	$(_end_touch)
rem_cyc_graph: $(REM_CYCLES_GRAPH)

GRAPH_STACK?=$(GRAPH_OUT_DIR2)/'.graphstack'
$(GRAPH_STACK):
	$(_start)
	mkdir -p $(GRAPH_OUT_DIR)/removed_cycles
	mkdir -p $(GRAPH_OUT_DIR)/removed_cycles/PCE_all_subs
	Rscript $(_r)/graph_long.r $(num_sub) $(FIG_OUT_DIR)/tablepath_ori_removed $(GRAPH_OUT_DIR) $(GRAPH_OUT_DIR)/PCE_all_subs $(_k) $(GRAPH_OUT_DIR)/med_cov_table
	mkdir -p $(GRAPH_OUT_DIR)/PCE_all_subs
	Rscript $(_r)/graph_long.r $(num_sub) $(FIG_OUT_DIR)/tablepath_ori $(GRAPH_OUT_DIR) $(GRAPH_OUT_DIR)/PCE_all_subs $(_k) $(GRAPH_OUT_DIR)/med_cov_table
graph_stack: $(GRAPH_STACK)

MED_COV?=$(GRAPH_OUT_DIR)/'.med_covdone'
$(MED_COV):
	python $(_py)/make_cov_table.py $(GRAPH_OUT_DIR) $(GRAPH_OUT_DIR)/med_cov_table
med_cov: $(MED_COV)

ALL_GRAPH?=$(GRAPH_OUT_DIR)/'.all_done'
$(ALL_GRAPH):
	Rscript $(_r)/all_pce_and_subs.r $(FIG_OUT_DIR)/tablepath_ori $(GRAPH_OUT_DIR) $(GRAPH_OUT_DIR)/med_cov_table
all_graph: $(ALL_GRAPH)
