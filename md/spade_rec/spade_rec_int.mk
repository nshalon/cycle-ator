#####################################################################################################
# register module
#####################################################################################################

units=sub_module_a.mk
$(call _register_module,spade_rec,$(units),,)

#####################################################################################################
# module parameters
#####################################################################################################

sample?=SRR059354

#read1=$(_read_dir)/$(sample)_1.fastq.gz
#read2=$(_read_dir)/$(sample)_2.fastq.gz
OUT_DIR=$(SAMPLE_DIR)/$(sample)
FIG_OUT_DIR=$(BASE_FDIR)/nobub_$(sample)
GRAPH_OUT_DIR=$(BASE_FDIR)/graphs/nobub_$(sample)
read1?=$(SAMPLE_DIR)/$(sample)/R1.fastq
read2?=$(SAMPLE_DIR)/$(sample)/R2.fastq

sample2?='dontrun'

ifdef ncbi
read1=/relman04/data/public/Mehta_2018/$(sample)_1.fastq
read2=/relman04/data/public/Mehta_2018/$(sample)_2.fastq
2_read1=/relman04/data/public/Mehta_2018/$(sample2)_1.fastq
2_read2=/relman04/data/public/Mehta_2018/$(sample2)_2.fastq
endif

FIG_OUT_DIR2=$(BASE_FDIR)/$(sample)/$(sample2)
GRAPH_OUT_DIR2=$(BASE_FDIR)/graphs/$(sample)/$(sample2)
