**SNP_calling_ASE.md** contains steps to call informative variants from RNAseq reads for allele-specific expression analyses (informative to MANA vs SARA ancestry). These steps generate the input for count_ASE_reads.md 

**count_ASE_reads.md** contains steps to align F1 reads with WASP to bin them into MANA and SARA groups, and then count the resulting reads using HTSeq. These steps generate the input for ASE analyses using DESseq2 in R.

**read_processing_counting.md** contains steps to align all pureline reads and count using HTSeq. These steps generate the input for differential expression analyses using DESseq2 in R.
