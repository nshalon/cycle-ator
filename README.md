# cCycle-ator

Cycleator is a tool for *de novo* recognition of circular mobile elements such as plasmids and phages in microbial metagenomic samples. It's designed with the purpose of avoiding false positives from chromosomal cycles that can arise in assembly graphs.

To access this pipeline, type

`$ git clone https://github.com/nshalon/cycleator.git`

This pipeline uses common bioinformatic software. Please edit the following variables in `global.mk` so that they lead to local paths for the following tools:

- \_mega: megahit for assembly  
- \_nucmer: nucmer for alignment (usually for simulated data)
- \_bwa: bwa for mapping reads

Other configurable variables:

- \_k: edit the kmer size used for assembly
- \_read\_dir: directory to find paired reads in fastq format called R1.fastq and R2.fastq (names can be edited in cyc\_find interface)

To run a sample, type

`$ make cyc_find sample=SAMPLE_NAME`


 
