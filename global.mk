###################################################
# TOOL AND SCRIPT  DIRECTORIES, EDIT AS NEEDED
###################################################

_recycle_dir?=/home/nshalon/python/python3/bin
_bwa?=/relman02/users/eitany/download/bwa-0.7.12/bwa
_py?=md/Python
_pl?=md/pl
_r?=md/R
_mega?=/home/dethlefs/bin
_nucmer?=/home/nshalon/work/pipe/mummer-4.0.0beta2
varisum_dir?=/relman02/users/eitany/ms_root/makeshift-modules/vari//bin.koch.stanford.edu/
_read_dir?=$(OUT_DIR)
# DEFAULT KMER SIZE FOR ASSEMBLY
_k?=77
sim?=F # whether sample is simulated

# DB FOR ANNOTATION
UNIREF_DIAMOND_DB_DIR=/relman03/work/users/eitany/bcc/diamond_db
GO_BASE_DIR=/relman03/work/users/eitany/bcc/GO

BASE_OUTPUT_DIR=$(CURDIR)/out_base
BASE_FDIR=$(CURDIR)/out_base/figures

ASSEMBLY_DIR=$(OUT_DIR)
GENES_FASTA_INPUT=$(CYC_OUT_DIR)/cycles.fasta


# define output directory
SAMPLE_DIR?=$(BASE_OUTPUT_DIR)/out

