all: pair_reads

DIR_DONE1=$(OUT_DIR_GEN)/.base_dir
$(DIR_DONE1):
	mkdir -p $(OUT_DIR_EX)
dir: $(DIR_DONE1)

MF_DONE?=$(OUT_DIR_GEN)/'fasta'
$(MF_DONE): $(DIR_DONE1)
	Rscript $(_r)/R_call.r $(_r)/gen_assembly.r gen.$(sample) \
           ofn.fasta=$(OUT_DIR_EX)/mol.fasta \
           ofn.table=$(OUT_DIR_EX)/mol.tab
	rc=$?; if [ $rc -ne 0 ]; then exit 1; fi
makefasta: $(MF_DONE)

COORD_DONE?=$(OUT_DIR_GEN)/'coords'
$(COORD_DONE): $(MF_DONE)
	Rscript $(_r)/R_call.r $(_r)/gen_read_coords.r gen.read.coords \
           genome.table.ifn=$(OUT_DIR_EX)/mol.tab \
           insert=200 insert.sd=100 read.length=150 \
           ofn=$(OUT_DIR_EX)/mol.coords
	rc=$?; if [ $rc -ne 0 ]; then exit 1; fi	
gencoords: $(COORD_DONE)

READ_PAIR?=$(OUT_DIR_GEN)/'readpairs'
$(READ_PAIR): $(COORD_DONE)
	perl $(_pl)/generate_reads_pair.pl $(OUT_DIR_EX)/mol.coords $(OUT_DIR_EX)/mol.fasta $(OUT_DIR_EX)/R1.fastq $(OUT_DIR_EX)/R2.fastq
	rc=$?; if [ $rc -ne 0 ]; then exit 1; fi
pair_reads: $(READ_PAIR)

