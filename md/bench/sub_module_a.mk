all: pair_reads

DIR_DONE=$(OUT_DIR_GEN)
$(DIR_DONE):
	mkdir -p $(OUT_DIR_GEN)
	cp $(fasta) $(OUT_DIR_GEN)/mol.fasta
	sed -i '1c\>m1' $(OUT_DIR_GEN)/mol.fasta 
dir: $(DIR_DONE)

TAB_DONE?=$(OUT_DIR_GEN)/'tab'
$(TAB_DONE): $(DIR_DONE)
	python $(_py)/genTab.py $(fasta) $(OUT_DIR_GEN)/mol.tab
maketab: $(TAB_DONE)

COORD_DONE?=$(OUT_DIR_GEN)/'coords'
$(COORD_DONE): $(TAB_DONE)
	Rscript $(_r)/R_call.r $(_r)/gen_read_coords.r gen.read.coords \
           genome.table.ifn=$(OUT_DIR_GEN)/mol.tab \
           insert=200 insert.sd=100 read.length=150 \
           ofn=$(OUT_DIR_GEN)/mol.coords
	rc=$?; if [ $rc -ne 0 ]; then exit 1; fi	
gencoords: $(COORD_DONE)

READ_PAIR?=$(OUT_DIR_GEN)/'readpairs'
$(READ_PAIR): $(COORD_DONE)
	perl $(_pl)/generate_reads_pair.pl $(OUT_DIR_GEN)/mol.coords $(OUT_DIR_GEN)/mol.fasta $(OUT_DIR_GEN)/R1.fastq $(OUT_DIR_GEN)/R2.fastq
	rc=$?; if [ $rc -ne 0 ]; then exit 1; fi
pair_reads: $(READ_PAIR)

