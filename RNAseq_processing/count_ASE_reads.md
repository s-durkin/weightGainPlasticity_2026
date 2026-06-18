# Count F1 reads for allele-specific expression analyses: Splits reads into MANA and SARA groups

## 1. Remap F1 reads with STAR v.2.7.11b providing the informative SNPs via WASP
```STAR --runMode alignReads``` \
```--runThreadN 16``` \
```--genomeDir /expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/STAR_genomedirectory``` \
```--readFilesCommand zcat``` \
```--readFilesIn ${Sample}_R1_cleaned.fq.gz ${Sample}_R2_cleaned.fq.gz``` \
```--outSAMtype BAM SortedByCoordinate``` \
```--varVCFfile /expanse/lustre/projects/uoa104/sdurkin/mat_effects/snp_calling/cleanvariants/BZ_NY_FixDiff.F1het.mm85.PHASE.vcf``` \
```--waspOutputMode SAMtag``` \
```--outSAMattributes NH HI AS nM NM MD vA vG vW``` \
```--outFilterMultimapNmax 1``` \
```--outFilterMismatchNmax 3``` \
```--outFileNamePrefix /expanse/lustre/projects/uoa104/sdurkin/mat_effects/F1_remapped_reads/${Sample}_F1remap\```



## 2. Separate the aligned bams into population specific pools based on genotype in the phased VCF
### Grep for passing WASP filtering (vW tag) and Grep for the appropriate vA:B tag (the allele tag denoting which allele the read maps to). NY reads are vA 1, BZ reads are vA 2**

**Command for NY/SARA reads**\
A.```samtools view ${Sample}_F1remapAligned.sortedByCoord.out.bam | grep "vW:i:1" | grep "vA:B:c,1"  | awk -F' ' '{if($18 ~ /2/){}else{print}}' > ${Sample}_F1remapAligned.sortedByCoord.NYreads.txt``` \
B.```samtools view -H ${Sample}_F1remapAligned.sortedByCoord.out.bam > ${Sample}.header.txt``` \
C.```cat ${Sample}.header.txt ${Sample}_F1remapAligned.sortedByCoord.NYreads.txt > ${Sample}_F1remapAligned.sortedByCoord.NYreads.bam``` 
<br><br>
**Command for BZ/MANA reads**\
A. ```samtools view ${Sample}_F1remapAligned.sortedByCoord.out.bam | grep "vW:i:1" | grep "vA:B:c,2" | awk -F' ' '{if($18 ~ /1/){}else{print}}' > ${Sample}_F1remapAligned.sortedByCoord.BZreads.txt``` \
B. ```cat ${Sample}.header.txt ${Sample}_F1remapAligned.sortedByCoord.BZreads.txt > ${Sample}_F1remapAligned.sortedByCoord.BZreads.bam``` 



## 3. Sort and index read separate bam files with samtools
**Command for NY/SARA reads**\
A. ```samtools sort ${Sample}_F1remapAligned.sortedByCoord.NYreads.bam -o ${Sample}_F1remapAligned.NYreads.sorted.bam``` \
B. ```samtools index ${Sample}_F1remapAligned.NYreads.sorted.bam```
<br><br>
**Command for BZ/MANA reads**\
A. ```samtools sort ${Sample}_F1remapAligned.sortedByCoord.BZreads.bam -o ${Sample}_F1remapAligned.BZreads.sorted.bam``` \
B. ```samtools index ${Sample}_F1remapAligned.BZreads.sorted.bam```



## 4. Count population specific F1 bam files with HTSeq v.0.11.0
**Command for NY/SARA reads**\
```htseq-count -f bam --order=pos --stranded=no ${Sample}_F1remapAligned.NYreads.sorted.bam ```\
```/expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/Mus_musculus.GRCm39.110.gtf ```\
```> /expanse/lustre/projects/uoa104/sdurkin/mat_effects/F1_ASE_counts/${Sample}.NYreads.count```
<br><br>
**Command for BZ/MANA reads**\
```htseq-count -f bam --order=pos --stranded=no ${Sample}_F1remapAligned.BZreads.sorted.bam ```\
```/expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/Mus_musculus.GRCm39.110.gtf ```\
```> /expanse/lustre/projects/uoa104/sdurkin/mat_effects/F1_ASE_counts/${Sample}.BZreads.count```



## 5. Merge counts with merge_tables.py
```python merge_tables.py F1_ASE_samples.txt > F1_ASE_samples_counts.txt ```


<br><br>
This is the input for ASE DESeq analyses
