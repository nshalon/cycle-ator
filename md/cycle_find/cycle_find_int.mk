#####################################################################################################
# register module
#####################################################################################################

units=sub_module_real.mk
ifeq ($(sim), T)
	units=sub_module_sim.mk
endif

$(call _register_module,cycle_find,$(units),,)

#####################################################################################################
# module parameters
#####################################################################################################
CYC_FIND_OUT_DIR=$(BASE_OUTPUT_DIR)/cyc_find_out
CYC_FIND_OUT=$(CYC_FIND_OUT_DIR)/$(sample)
read1?=$(_read_dir)/R1.fastq
read2?=$(_read_dir)/R2.fastq
molfasta?=$(OUT_DIR)/mol.fasta
origfasta?=$(OUT_DIR)/ref.fa

