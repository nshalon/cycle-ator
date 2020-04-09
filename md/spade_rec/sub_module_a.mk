all: nucmer

DIR_DONE?=$(OUT_DIR)
$(DIR_DONE):
	mkdir -p $(OUT_DIR)
dir: $(DIR_DONE)

ASSEMBLE_DONE?=$(OUT_DIR)/.assemble_done
$(ASSEMBLE_DONE): $(DIR_DONE)
	$(_start)
	#/usr/bin/time -v -o $(OUT_DIR)/assembly_time $(_mega)/megahit -m 0.84 -o $(OUT_DIR)/megahit --min-contig-len 300 --k-min 27 --k-max 99 --k-step 10 --merge-level 20,0.95 -1 $(read1) -2 $(read2) -t 40
	/usr/bin/time -v -o $(OUT_DIR)/assembly_time $(_mega)/megahit -m 0.84 -o $(OUT_DIR)/megahit --min-contig-len 100 --k-min 27 --k-max 99 --k-step 10 --merge-level 1000,0.95 -1 $(read1) -2 $(read2) -t 40
	#$(_mega)/megahit -m 0.5 -o $(OUT_DIR)/megahit --min-contig-len 300 --k-list 147 --merge-level 20,0.95 -t 40 -1 $(read1) -2 $(read2)
	#python $(_spades) -k 55 --phred-offset 33 -1 $(read1) -2 $(read2) -o $(OUT_DIR)/assembly
	$(_end_touch)
assemble: $(ASSEMBLE_DONE)

FASTA_TO_FASTG_DONE?=$(OUT_DIR)/.fasta_to_fastg
$(FASTA_TO_FASTG_DONE): $(ASSEMBLE_DONE)
	$(_start)
	cp $(OUT_DIR)/megahit/intermediate_contigs/k$(_k).contigs.fa $(OUT_DIR)/megahit/final.fa
	$(_mega)/megahit_toolkit contig2fastg $(_k) $(OUT_DIR)/megahit/intermediate_contigs/k$(_k).contigs.fa > $(OUT_DIR)/k$(_k).fastg
	$(_end_touch)
fasta_to_fastg: $(FASTA_TO_FASTG_DONE)

BWA_DONE?=$(OUT_DIR)/.bwa_done
$(BWA_DONE): $(FASTA_TO_FASTG_DONE)
	$(_start)
	/usr/bin/time -v -o $(OUT_DIR)/bwa_index_time $(_bwa) index $(OUT_DIR)/megahit/final.fa
	$(_end_touch)
bwa: $(BWA_DONE)

SAM_ONE_DONE?=$(OUT_DIR)/.sam_one_done
$(SAM_ONE_DONE): $(BWA_DONE)
	$(_start)
	/usr/bin/time -v -o $(OUT_DIR)/bwa_mem_time $(_bwa) mem -t 40 $(OUT_DIR)/megahit/final.fa $(read1) $(read2) | samtools view -buS - > $(OUT_DIR)/reads_pe.bam
	$(_end_touch)
sam_one: $(SAM_ONE_DONE)

SAM_TWO_DONE?=$(OUT_DIR)/.sam_two_done
$(SAM_TWO_DONE): $(SAM_ONE_DONE)
	$(_start)
	samtools view --threads 40 -bF 0x0800 $(OUT_DIR)/reads_pe.bam > $(OUT_DIR)/reads_pe_primary.bam
	samtools sort --threads 40 $(OUT_DIR)/reads_pe_primary.bam -o $(OUT_DIR)/Sorted.BAM
	samtools index $(OUT_DIR)/Sorted.BAM
	$(_end_touch)
sam_two: $(SAM_TWO_DONE)

RECY_DONE?=$(OUT_DIR)/.recy_done
$(RECY_DONE): $(SAM_TWO_DONE)
	$(_start)
	mkdir -p $(FIG_OUT_DIR)
	mkdir -p $(OUT_DIR)/unfil_recycle
	python $(_recycle_dir)/edit_recycle.py -g $(OUT_DIR)/k$(_k).fastg -k $(_k) -m 100 -b $(OUT_DIR)/Sorted.BAM -i False -o $(OUT_DIR)/unfil_recycle
	$(_end_touch)
recycler: $(RECY_DONE)

RENAME_MEGA_FASTA?=$(FIG_OUT_DIR)/.rename_mega_done
$(RENAME_MEGA_FASTA): $(RECY_DONE)
	$(_start)
	python $(_py)/rename_megahit_fasta.py $(OUT_DIR)/megahit/final.fa $(OUT_DIR)/megahit/rename_final.contigs.fa $(_k)
	$(_end_touch)
rename_mega: $(RENAME_MEGA_FASTA)


#REMOVE_LONG_PCE?=$(FIG_OUT_DIR)/'removed_long_pce'
#$(REMOVE_LONG_PCE): $(RENAME_MEGA_FASTA)
#	$(_start)
#	python $(_py)/remove_large_CE.py $(OUT_DIR)/unfil_recycle/k$(_k).cycs.paths_w_cov.txt $(OUT_DIR)/unfil_recycle/edit_k$(_k).cycs.paths_w_cov.txt
#	$(_end_touch)
#remove_long: $(REMOVE_LONG_PCE)

RENAME_DONE?=$(FIG_OUT_DIR)/.rename_done
$(RENAME_DONE): $(RENAME_MEGA_FASTA)
	$(_start)
	mkdir -p $(FIG_OUT_DIR)
	python $(_py)/rename_path.py $(OUT_DIR)/unfil_recycle/k$(_k).cycs.paths_w_cov.txt $(OUT_DIR)/megahit/rename_final.contigs.fa $(FIG_OUT_DIR)/renamefasta $(FIG_OUT_DIR)/tablepath $(_k)
	$(_end_touch)
rename: $(RENAME_DONE)

ORI_DONE?=$(FIG_OUT_DIR)/.ori_done
$(ORI_DONE): $(RENAME_DONE)
	$(_start)
	python $(_py)/orientation_v2.py $(FIG_OUT_DIR)/renamefasta $(OUT_DIR)/unfil_recycle/k$(_k).cycs.fasta  $(FIG_OUT_DIR)/tablepath $(FIG_OUT_DIR)/tablepath_ori $(_k)
	$(_end_touch)
orientation: $(ORI_DONE)

TEMP_SUM?=$(GRAPH_OUT_DIR)/.tempsum
$(TEMP_SUM): $(ORI_DONE)
	$(_start)
	mkdir -p $(GRAPH_OUT_DIR)
	python $(_py)/tempor_sum.py $(FIG_OUT_DIR)/tablepath_ori $(GRAPH_OUT_DIR)/summary $(_k)
	$(_end_touch)
temp_sum: $(TEMP_SUM)

MAP_DONE?=$(FIG_OUT_DIR)/.map_done
$(MAP_DONE): $(TEMP_SUM)
	$(_start)
	/usr/bin/time -v -o $(OUT_DIR)/bwa_index2_time $(_bwa) index $(FIG_OUT_DIR)/renamefasta
	gunzip -c -f $(read1) > $(FIG_OUT_DIR)/R1.fastq
	gunzip -c -f $(read2) > $(FIG_OUT_DIR)/R2.fastq
	perl $(_pl)/trim_fastq.pl $(FIG_OUT_DIR)/R1.fastq 0 40 $(FIG_OUT_DIR)/trimR1.fastq
	perl $(_pl)/trim_fastq.pl $(FIG_OUT_DIR)/R2.fastq 0 40 $(FIG_OUT_DIR)/trimR2.fastq
	/usr/bin/time -v -o $(FIG_OUT_DIR)/map_one_side_time $(_bwa) mem $(FIG_OUT_DIR)/renamefasta $(FIG_OUT_DIR)/trimR1.fastq > $(FIG_OUT_DIR)/r1-aln.sam
	$(_bwa) mem $(FIG_OUT_DIR)/renamefasta $(FIG_OUT_DIR)/trimR2.fastq > $(FIG_OUT_DIR)/r2-aln.sam
	perl $(_pl)/parse_bwa_sam.pl $(FIG_OUT_DIR)/r1-aln.sam $(FIG_OUT_DIR)/R1table $(FIG_OUT_DIR)/R1table_stats
	perl $(_pl)/parse_bwa_sam.pl $(FIG_OUT_DIR)/r2-aln.sam $(FIG_OUT_DIR)/R2table $(FIG_OUT_DIR)/R2table_stats
	/usr/bin/time -v -o $(FIG_OUT_DIR)/filter_one_side_time perl $(_pl)/filter_map.pl $(FIG_OUT_DIR)/R1table 30 40 2 $(FIG_OUT_DIR)/filter_R1table $(FIG_OUT_DIR)/filter_R1table_stats
	perl $(_pl)/filter_map.pl $(FIG_OUT_DIR)/R2table 30 40 2 $(FIG_OUT_DIR)/filter_R2table $(FIG_OUT_DIR)/filter_R2table_stats
	$(_end_touch)
map: $(MAP_DONE)

PAIR_READS?=$(FIG_OUT_DIR)/'.pair_done'
$(PAIR_READS): $(MAP_DONE)
	$(_start)	
	/usr/bin/time -v -o $(FIG_OUT_DIR)/pair_reads_time perl $(_pl)/pair_reads.pl $(FIG_OUT_DIR)/filter_R1table $(FIG_OUT_DIR)/filter_R2table $(FIG_OUT_DIR)/paired_table_pre $(FIG_OUT_DIR)/paired_table_stats
	$(_end_touch)
pair_reads: $(PAIR_READS)

REM_NCYC?=$(FIG_OUT_DIR)/'.rem_done'
$(REM_NCYC): $(PAIR_READS)
	$(_start)
	python $(_py)/rem_noncyc_reads.py $(FIG_OUT_DIR)/paired_table_pre $(FIG_OUT_DIR)/paired_table $(FIG_OUT_DIR)/paired_table_remove_stats
	$(_end_touch)
remove_reads: $(REM_NCYC)

#ORI_DONE?=$(FIG_OUT_DIR)/'.ori_done'
#$(ORI_DONE): $(REM_NCYC)
#	$(_start)
#	python $(_py)/orientation.py $(FIG_OUT_DIR)/renamefasta $(OUT_DIR)/unfil_recycle/k$(_k).cycs.fasta $(FIG_OUT_DIR)/paired_table $(FIG_OUT_DIR)/tablepath $(FIG_OUT_DIR)/tablepath_ori $(_k)
#	$(_end_touch)
#orientation: $(ORI_DONE)

CLASS_DONE?=$(FIG_OUT_DIR)/'class_done'
$(CLASS_DONE): $(REM_NCYC)
	$(_start)
	mkdir -p $(FIG_OUT_DIR)/pairedtables
	/usr/bin/time -v -o $(FIG_OUT_DIR)/classify_reads_time python $(_py)/map_contigs.py $(FIG_OUT_DIR)/paired_table $(FIG_OUT_DIR)/tablepath_ori 2000 $(FIG_OUT_DIR)/classified_reads $(FIG_OUT_DIR)/pairedtables $(_k) $(FIG_OUT_DIR)/class_stats
	$(_end_touch)
class: $(CLASS_DONE)

BUCK_DONE?=$(FIG_OUT_DIR)/'bucket_done'
$(BUCK_DONE): $(CLASS_DONE)
	$(_start)
	python $(_py)/table_bucket.py $(FIG_OUT_DIR)/classified_reads 100 $(FIG_OUT_DIR)/bucket_classified_reads_100
	python $(_py)/table_bucket.py $(FIG_OUT_DIR)/classified_reads 10 $(FIG_OUT_DIR)/bucket_classified_reads_10
	$(_end_touch)
bucket: $(BUCK_DONE)

PLOT_DONE?=$(FIG_OUT_DIR)/'plot_done'
$(PLOT_DONE): $(BUCK_DONE)
	$(_start)
	mkdir -p $(GRAPH_OUT_DIR)
	mkdir -p $(GRAPH_OUT_DIR)/bucket100
	mkdir -p $(GRAPH_OUT_DIR)/bucket10
	Rscript $(_r)/plot_table.r $(FIG_OUT_DIR)/bucket_classified_reads_100 $(FIG_OUT_DIR)/tablepath_ori 100 $(GRAPH_OUT_DIR)/bucket100
	Rscript $(_r)/plot_table.r $(FIG_OUT_DIR)/bucket_classified_reads_10 $(FIG_OUT_DIR)/tablepath_ori 10 $(GRAPH_OUT_DIR)/bucket10
	$(_end_touch)
graph: $(PLOT_DONE)

COV_DONE?=$(GRAPH_OUT_DIR)/'cov_done'
$(COV_DONE): $(PLOT_DONE)
	$(_start)
	mkdir -p $(GRAPH_OUT_DIR)/cov
	/usr/bin/time -v -o $(FIG_OUT_DIR)/cov_time python $(_py)/xcov.py 1 $(FIG_OUT_DIR)/paired_table $(FIG_OUT_DIR)/tablepath_ori $(GRAPH_OUT_DIR)/cov $(GRAPH_OUT_DIR)/cov_table $(GRAPH_OUT_DIR)/stats $(_k)
	$(_end_touch)
coverage: $(COV_DONE)

GRAPH_COV?=$(GRAPH_OUT_DIR)/'graph_done'
$(GRAPH_COV): $(COV_DONE)
	$(_start)
	Rscript $(_r)/graph_cov.R $(GRAPH_OUT_DIR)/cov_table 1 $(FIG_OUT_DIR)/tablepath_ori $(GRAPH_OUT_DIR)/cov $(GRAPH_OUT_DIR)/stats $(_k)
	$(_end_touch)
graph_cov: $(GRAPH_COV)
	
SUM_DONE?=$(FIG_OUT_DIR)/'summary_done'
$(SUM_DONE): $(GRAPH_COV)
	$(_start)
	rm $(GRAPH_OUT_DIR)/summary
	python $(_py)/sample_sum.py $(OUT_DIR)/unfil_recycle/k$(_k)_CV_profile.txt $(FIG_OUT_DIR)/tablepath_ori $(FIG_OUT_DIR)/bucket_classified_reads_100 $(GRAPH_OUT_DIR)/summary $(_k)
	$(_end_touch)
summary: $(SUM_DONE)

NUCMER=$(FIG_OUT_DIR)/nucmer.done
$(NUCMER): #$(SUM_DONE)
	$(_nucmer)/nucmer --prefix=REF_REF $(fasta) $(fasta)
	mv /home/nshalon/work/pipe/REF_REF.delta $(FIG_OUT_DIR)
	$(_nucmer)/nucmer --prefix=REF_CONTIG $(fasta) $(OUT_DIR)/megahit/final.contigs.fa
	mv /home/nshalon/work/pipe/REF_CONTIG.delta $(FIG_OUT_DIR)
	$(_nucmer)/nucmer --prefix=CONTIG_CONTIG $(OUT_DIR)/megahit/final.contigs.fa $(OUT_DIR)/megahit/final.contigs.fa
	mv /home/nshalon/work/pipe/CONTIG_CONTIG.delta $(FIG_OUT_DIR)
	$(_nucmer)/mummerplot --prefix=REF_REF --png $(FIG_OUT_DIR)/REF_REF.delta
	$(_nucmer)/mummerplot --prefix=REF_CONTIG --png $(FIG_OUT_DIR)/REF_CONTIG.delta
	$(_nucmer)/mummerplot --prefix=CONTIG_CONTIG --png $(FIG_OUT_DIR)/CONTIG_CONTIG.delta
	mv REF_* $(GRAPH_OUT_DIR)
	mv CONTIG_* $(GRAPH_OUT_DIR)
nucmer: $(NUCMER)
