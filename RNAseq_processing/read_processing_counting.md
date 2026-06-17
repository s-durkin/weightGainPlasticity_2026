## RNAseq read initial processing and counting

# 1. Quality Filter and Adapter Trimming Using FastP v.0.19.6

```/expanse/lustre/projects/uoa104/sdurkin/fastp/fastp \```\
        ```-i $seq\_R1_001.fastq.gz \```\
        ```-I $seq\_R2_001.fastq.gz \```\
        ```-o $seq\_R1_cleaned.fq.gz \```\
        ```-O $seq\_R2_cleaned.fq.gz \```\
        ```-n 5 \```\
        ```-q 15 \```\
        ```-u 30 \```\
        ```--detect_adapter_for_pe \```\
        ```--cut_window_size=4 \```\
        ```--cut_mean_quality=15 \```\
        ```--length_required=25 \```\
        ```-j $seq\_report.json \```\
        ```-h $seq\_report.html \```\
        ```-w 2;
```mv -v "${seq}_R1_cleaned.fq.gz" "../clean_reads/${seq}_R1_cleaned.fq.gz"```\
```mv -v "${seq}_R2_cleaned.fq.gz" "../clean_reads/${seq}_R2_cleaned.fq.gz"```


# 2. Align clean reads with STAR v.2.7.11b

**FIRST: Create genome index - GRCm39 build** 
```STAR --runThreadN 16 \```\
```--runMode genomeGenerate \```\
```--limitGenomeGenerateRAM 33524399488 \```\
```--genomeDir /expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/STAR_genomedirectory \```\
```--genomeFastaFiles /expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa \```\
```--sjdbGTFfile /expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/Mus_musculus.GRCm39.110.gtf \```\
```--sjdbOverhang 149```

**SECOND: Run STAR**
```STAR --runMode alignReads \```\
```--runThreadN 16 \```\
```--genomeDir /expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/STAR_genomedirectory \```\
```--readFilesCommand zcat \```\
```--readFilesIn ${Sample}_R1_cleaned.fq.gz ${Sample}_R2_cleaned.fq.gz \```\
```--outSAMtype BAM SortedByCoordinate \```\
```--outFilterMultimapNmax 1 \```\
```--outFilterMismatchNmax 3 \```\
```--outFileNamePrefix /expanse/lustre/projects/uoa104/sdurkin/mat_effects/mapped_reads/${Sample}```


# 3. Count aligned reads with HTSeq v.0.11.0

```htseq-count -f bam \```\
```--order=pos \```\
```--stranded=no \```\
```${Sample}Aligned.sortedByCoord.out.bam \```\
```/expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/Mus_musculus.GRCm39.110.gtf \```\
```> /expanse/lustre/projects/uoa104/sdurkin/mat_effects/countdata/${line}.liver.count```


# 4. Merge counts with merge_tables.py
```python merge_tables.py all_samples.txt > all_samples_counts.txt```


