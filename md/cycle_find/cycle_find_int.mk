#####################################################################################################
# register module
#####################################################################################################

units=sub_module_a.mk
$(call _register_module,cycle_find,$(units),,)

#####################################################################################################
# module parameters
#####################################################################################################

CYC_FIND_OUT_DIR=$(BASE_OUTPUT_DIR)/cycle_find
CYC_FIND_OUT=$(CYC_FIND_OUT_DIR)/$(sample)
read1?=$(OUT_DIR)/R1.fastq
read2?=$(OUT_DIR)/R2.fastq
molfasta?=$(OUT_DIR)/mol.fasta
origfasta?=$(OUT_DIR)/ref.fa

