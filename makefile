include $(MAKESHIFT_ROOT)/makeshift-core/makeshift.mk

#####################################################################################################
# config file
#####################################################################################################

c?=config/default.cfg
include $(c)
$(call _set_user_title,config: $(c))

#####################################################################################################
# modules
#####################################################################################################

# global parameters
$(call _module,global.mk)
$(call _module_root,metagenomics/genes)

# include module
$(call _module,md/gen_ex/gen_ex_int.mk)
$(call _module,md/spade_rec/spade_rec_int.mk)
$(call _module,md/annotation/annot_int.mk)
$(call _module,md/map_long/map_long_int.mk)
$(call _module,md/bench/bench_int.mk)
$(call _module,md/cycle_find/cycle_find_int.mk)

# default module
m?=gen_ex
$(call _active_module,$(m))

#####################################################################################################
# main pipeline rules
#####################################################################################################
bench:
	@$(MAKE) m=bench pair_reads
make_ex:
	@$(MAKE) m=gen_ex pair_reads 
spades:
	@$(MAKE) m=spade_rec all
pgenes:
	@$(MAKE) m=genes prodigal genes_uniref genes_GO
	@$(MAKE) m=annotation tables
gene_table:
	@$(MAKE) m=annotation cov_genes
mlong:
	@$(MAKE) m=map_long all
cyc_find:
	@$(MAKE) m=cycle_find all
