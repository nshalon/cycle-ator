MAKE_TABLE?=$(GRAPH_OUT_DIR)/'.genes_done'
$(MAKE_TABLE):
	python $(_py)/make_gene_tables.py $(prodigal_table) $(uniref_table) /relman03/work/users/nshalon/pipe/out_base/figures/graphs/$(sample)/summary /relman03/work/users/nshalon/pipe/out_base/figures/graphs/$(sample)/gene_table.csv 
tables: $(MAKE_TABLE)

HMM_TABLE?=$(GRAPH_OUT_DIR)/'.hmm_done'
$(HMM_TABLE):
	python $(_py)/make_hmm_table.py $(prodigal_table) /home/nshalon/work/pipe/hmm_profiles/out $(GRAPH_OUT_DIR)/summary $(GRAPH_OUT_DIR)/hmm.csv $(GRAPH_OUT_DIR)/hmm_matrix.csv
hmm: $(HMM_TABLE)

PARSE_PFAM?=$(GRAPH_OUT_DIR)/'.parse_pfam'
$(PARSE_PFAM):
	python $(_py)/parse_and_classify_pfams.py $(prodigal_table) $(pfam) $(GRAPH_OUT_DIR)/summary $(GRAPH_OUT_DIR)/pfam_gene_table.csv $(GRAPH_OUT_DIR)/pfam_matrix.csv $(GRAPH_OUT_DIR)/pfam_gene_class.csv
parse_pfam: $(PARSE_PFAM)

MAKE_COV_W_GENES?=$(GRAPH_OUT_DIR)/'.genes'
$(MAKE_COV_W_GENES):
	mkdir -p $(GRAPH_OUT_DIR)/cov_with_pfam_genes
	mkdir -p $(GRAPH_OUT_DIR)/cov_with_genes_gc
	Rscript $(_r)/cov_with_gene_table.r $(FIG_OUT_DIR)/tablepath_ori $(GRAPH_OUT_DIR) $(GRAPH_OUT_DIR)/cov_with_genes_gc $(_k) $(GRAPH_OUT_DIR)/med_cov_table $(GRAPH_OUT_DIR)/gene_table.csv /home/nshalon/work/pipe/pfam_gene_class_temp.csv /home/nshalon/work/pipe/gc_file 
cov_genes: $(MAKE_COV_W_GENES)

	
