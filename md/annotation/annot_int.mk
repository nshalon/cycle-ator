units=sub_module_a.mk
$(call _register_module,annotation,$(units),,)
        
#####################################################################################################
# module parameters
# #####################################################################################################
#     
sample?=SRR059354

GENE_DIR=/relman03/work/users/nshalon/pipe/out_base/figures/graphs/$(sample)/genes
#
pfam?=/home/nshalon/work/pipe/hmm_profiles/db/out/Pfam_out_eval_over_5
prodigal_table=$(GENE_TABLE)
uniref_table=$(SAMPLE_DIR)/$(sample)/uniref/2019_03/table_uniq_taxa
