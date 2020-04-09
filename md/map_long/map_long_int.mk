units=sub_module_a.mk
$(call _register_module,map_long,$(units),,)


2_read1=/relman04/data/public/Mehta_2018/$(sample2)_1.fastq
2_read2=/relman04/data/public/Mehta_2018/$(sample2)_2.fastq

ifdef agg10
2_read1=/relman04/data/public/Mehta_2018/agg_10/$(sample2)_1.fastq
2_read2=/relman04/data/public/Mehta_2018/agg_10/$(sample2)_2.fastq
endif

FIG_OUT_DIR2=$(BASE_FDIR)/$(sample)/$(sample2)
GRAPH_OUT_DIR2=$(BASE_FDIR)/graphs/$(sample)/$(sample2)

_k?=147

FIG_OUT_DIR=$(BASE_FDIR)/$(sample)
GRAPH_OUT_DIR=$(BASE_FDIR)/graphs/$(sample)

num_sub?=4

sample2?='dontrun'
