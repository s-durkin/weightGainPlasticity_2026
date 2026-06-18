# SNP calling for allele-specific expression analyses


## 1. Mark Duplicates with picard version 2.20.8
**Done for pureline samples and F1 samples**
```java -jar $PICARDHOME/bin/picard.jar MarkDuplicates```\
```I=/expanse/lustre/projects/uoa104/sdurkin/mat_effects/mapped_reads/${Sample}Aligned.sortedByCoord.out.bam``` \
```O=${Sample}.sortedByCoord.markdups.bam``` \
```M=${Sample}metrics.txt ```



## 2. Split N Cigar reads with GATK version 4.2.2.0
**Done for pureline samples and F1 samples**
```gatk SplitNCigarReads``` \
```-R /expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa``` \
```-I ${Sample}.sortedByCoord.markdups.bam``` \
```-O ${Sample}.sortedByCoord.markdups.split.bam```



## 3. Add read groups with with picard version 2.20.8
**Done for pureline samples and F1 samples**
```java -jar $PICARDHOME/bin/picard.jar AddOrReplaceReadGroups``` \
```I=${line1}.sortedByCoord.markdups.split.bam``` \
```O=${line1}.sortedByCoord.markdups.RG.split.bam``` \
```RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${line2}```




## 4. Index with samtools
**Done for pureline samples and F1 samples**
```samtools index ${FILE}```



## 5. Call haplotypes with GATK HaplotypeCaller, version 4.2.2.0
**Done for pureline samples and F1 samples**
```gatk HaplotypeCaller -R /expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa``` \
```-I ${liver_IND}.sortedByCoord.markdups.RG.split.bam``` \
```-ERC GVCF --dont-use-soft-clipped-bases -stand-call-conf 20``` \
```-O /expanse/lustre/projects/uoa104/sdurkin/mat_effects/snp_calling/rawvariants/${sample_name}_rawvariants_combo.g.vcf.gz```



## 6. Combine GVCFs with GATK , version 4.2.2.0
**Command for pureline samples**
```gatk CombineGVCFs -R /expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa``` \
```--variant BB-18_rawvariants_combo.g.vcf.gz``` \
```--variant BB-19_rawvariants_combo.g.vcf.gz``` \
```--variant BB-20_rawvariants_combo.g.vcf.gz``` \
```--variant BB-21_rawvariants_combo.g.vcf.gz``` \
... \
```--variant NN-77_rawvariants_combo.g.vcf.gz``` \
```--variant NN-78_rawvariants_combo.g.vcf.gz``` \
```--variant NN-80_rawvariants_combo.g.vcf.gz``` \
```--variant NN-86_rawvariants_combo.g.vcf.gz``` \
```--O pureline.merged.rawvariants.g.vcf```
<br><br>
**Command for F1 samples**
```gatk CombineGVCFs -R /expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa``` \
```--variant BN-10_rawvariants_combo.g.vcf.gz``` \
```--variant BN-12_rawvariants_combo.g.vcf.gz``` \
```--variant BN-14_rawvariants_combo.g.vcf.gz``` \
```--variant BN-1_rawvariants_combo.g.vcf.gz``` \
... \
```--variant NB-3_rawvariants_combo.g.vcf.gz``` \
```--variant NB-6_rawvariants_combo.g.vcf.gz``` \
```--variant NB-7_rawvariants_combo.g.vcf.gz``` \
```--variant NB-8_rawvariants_combo.g.vcf.gz``` \
```--O F1s.merged.rawvariants.g.vcf```


## 7. Genotype GVCFs with GATK , version 4.2.2.0
gatk GenotypeGVCFs \
-R /expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa \
--variant pureline.merged.rawvariants.g.vcf \
--O pureline.merged.rawvariants.jointGenotyped.vcf.gz
Done for pureline and F1 separately
Steps 7-9 dones in one script: Geno_SNPfilter_pureline.sh and Geno_SNPfilter_F1s.sh


## 8. Select SNPs with GATK , version 4.2.2.0
gatk SelectVariants \
-R /expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa \
--variant pureline.merged.rawvariants.jointGenotyped.vcf.gz \
--select-type-to-include SNP -O pureline.merged.rawvariants.jointGenotyped.SNPs.vcf.gz
Done for pureline and F1 separately


## 9. Filter SNPs for high quality with GATK , version 4.2.2.0
gatk VariantFiltration \
-R /expanse/lustre/projects/uoa104/sdurkin/mus_resources/m39_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa \
--variant pureline.merged.rawvariants.jointGenotyped.SNPs.vcf.gz \
--filter-name "QD_2" --filter-expression "QD < 2.0" \
--filter-name "FS_60" --filter-expression "FS > 60.0" \
--filter-name "MQ_40" --filter-expression "MQ < 40.0" \
--filter-name "MQRankSum_-12.5" --filter-expression "MQRankSum < -12.5" \
--filter-name "ReadPosRankSum_-8" --filter-expression "ReadPosRankSum < -8.0" \
-O /expanse/lustre/projects/uoa104/sdurkin/mat_effects/snp_calling/cleanvariants/pureline.merged.rawvariants.jointGenotyped.SNPs.filter.vcf.gz
Done for pureline and F1 separately


## 10. Select only SNPs with PASS filter status using bcftools
bcftools view -i 'FILTER="PASS"' pureline.merged.rawvariants.jointGenotyped.SNPs.filter.vcf.gz \
-O z -o pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.vcf.gz
bcftools view -i 'FILTER="PASS"' F1s.merged.rawvariants.jointGenotyped.SNPs.filter.vcf.gz \
-O z -o F1s.merged.rawvariants.jointGenotyped.SNPs.PASS.vcf.gz
filter_PASS.sh


## 11. PURELINE: Create separate BZ and NY VCFs for pureline samples and filter for 10/12 individuals having data
bcftools view -s BB-18,BB-19,BB-20,BB-21,BB-75,BB-77,BB-78,BB-79,BB-81,BB-82,BB-83,BB-84 pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.vcf.gz -Oz -o BZ_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.vcf.gz
tabix BZ_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.vcf.gz
vcftools --gzvcf BZ_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.vcf.gz --max-missing 0.85 --recode --recode-INFO-all --stdout | bgzip -c > BZ_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.vcf.gz
After filtering, kept 480396 out of a possible 3162791 Sites

bcftools view -s NN-10,NN-12,NN-58,NN-59,NN-61,NN-62,NN-66,NN-67,NN-77,NN-78,NN-80,NN-86 pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.vcf.gz -Oz -o NY_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.vcf.gz
tabix NY_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.vcf.gz
vcftools --gzvcf NY_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.vcf.gz --max-missing 0.85 --recode --recode-INFO-all --stdout | bgzip -c > NY_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.vcf.gz
After filtering, kept 458535 out of a possible 3162791 Sites


## 12. PURELINE: Create separate fixed and different VCF files using filtering for hom and ref reciprocally 
filter the line specific VCFs for fixed hom ref, and fixed hom alt SNPs (two files per line)
#homozygous alt - tabix after
vcftools --gzvcf BZ_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.vcf.gz --non-ref-af 1 --recode --recode-INFO-all --stdout | bgzip -c > BZ_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.HomALT.vcf.gz
vcftools --gzvcf NY_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.vcf.gz --non-ref-af 1 --recode --recode-INFO-all --stdout | bgzip -c > NY_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.HomALT.vcf.gz

#homozygous ref - tabix after
vcftools --gzvcf BZ_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.vcf.gz --max-non-ref-af 0 --recode --recode-INFO-all --stdout | bgzip -c > BZ_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.HomREF.vcf.gz
vcftools --gzvcf NY_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.vcf.gz --max-non-ref-af 0 --recode --recode-INFO-all --stdout | bgzip -c > NY_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.HomREF.vcf.gz

Merge to create a fixed and different SNP file - tabix after
vcf-merge BZ_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.HomALT.vcf.gz NY_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.HomREF.vcf.gz | bgzip -c > BZ_NY_alt.ref.vcf.gz
vcftools --gzvcf BZ_NY_alt.ref.vcf.gz --max-missing 0.85 --recode --recode-INFO-all --stdout | bgzip -c > BZ_NY_alt.ref.mm85.vcf.gz
kept 89478 out of a possible 354301 Sites

vcf-merge BZ_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.HomREF.vcf.gz NY_pureline.merged.rawvariants.jointGenotyped.SNPs.PASS.mm85.HomALT.vcf.gz | bgzip -c > BZ_NY_ref.alt.vcf.gz
vcftools --gzvcf BZ_NY_ref.alt.vcf.gz --max-missing 0.85 --recode --recode-INFO-all --stdout | bgzip -c > BZ_NY_ref.alt.mm85.vcf.gz
kept 66974 out of a possible 312600 Sites

Final merge - just combine them together
bcftools concat --allow-overlaps BZ_NY_alt.ref.mm85.vcf.gz BZ_NY_ref.alt.mm85.vcf.gz -Oz -o BZ_NY_fixDiff.mm85.vcf.gz


## 13. F1s: Select heterzygous SNPs using GATK VariantFiltration isHet and filter for 20/24 individuals having data
gatk VariantFiltration \
-V F1s.merged.rawvariants.jointGenotyped.SNPs.PASS.vcf.gz \
-O F1s.merged.rawvariants.jointGenotyped.SNPs.isHet.vcf.gz \
--genotype-filter-name "isHetFilter" --genotype-filter-expression "isHet == 1" \
--genotype-filter-name "isHomRefFilter" --genotype-filter-expression "isHomRef == 1" \
--genotype-filter-name "isHomAltFilter" --genotype-filter-expression "isHomVar == 1"
Script is F1s_isHetFilter.sh

Select just the "isHetFilter" snps
bcftools filter -e 'FORMAT/FT="isHomRefFilter"' F1s.merged.rawvariants.jointGenotyped.SNPs.isHet.vcf.gz -O z -o F1s.merged.rawvariants.jointGenotyped.SNPs.noHomRef.vcf.gz
bcftools filter -e 'FORMAT/FT="isHomAltFilter"' F1s.merged.rawvariants.jointGenotyped.SNPs.noHomRef.vcf.gz -O z -o F1s.merged.rawvariants.jointGenotyped.SNPs.noHomRefAlt.vcf.gz
bcftools filter -i 'FORMAT/FT="isHetFilter"' F1s.merged.rawvariants.jointGenotyped.SNPs.noHomRefAlt.vcf.gz -O z -o F1s.merged.rawvariants.jointGenotyped.SNPs.isHetPASS.vcf.gz This step should be redundant

Filter for missingness
vcftools --gzvcf F1s.merged.rawvariants.jointGenotyped.SNPs.isHetPASS.vcf.gz --max-missing 0.85 --recode --recode-INFO-all --stdout | bgzip -c > F1hetCalls.mm85.vcf.gz
After filtering, kept 67039 out of a possible 72112 Sites


## 14. Merge F1 het VCF and fixed and different pureline VCF and filter for 80% having data
vcf-merge BZ_NY_fixDiff.mm85.vcf.gz F1hetCalls.mm85.vcf.gz | bgzip -c > BZ_NY_FixDiff.F1het.vcf.gz
vcftools --gzvcf BZ_NY_FixDiff.F1het.vcf.gz --max-missing 0.85 --recode --recode-INFO-all --stdout | bgzip -c > BZ_NY_FixDiff.F1het.mm85.vcf.gz
After filtering, kept 63280 out of a possible 160211 Sites
was more like 25000-30000 for my first F1 project, so much better!



## 15. Phase Fix different/F1 het VCF to have  0|1 or 1|0 corresponding to whether the alt is in BZ or NY, respectively 
This is what we use for WASP
gunzip -c BZ_NY_FixDiff.F1het.mm85.vcf.gz | awk -F' ' '{if($1 ~ /^#/){}else{print}}' | awk -F' ' '{if( ($10 ~/^1\/1/ || $10 ~/^1\|1/) || ($11 ~/^1\/1/ || $11 ~/^1\|1/) || ($12 ~/^1\/1/ || $12 ~/^1\|1/) || ($13 ~/^1\/1/ || $13 ~/^1\|1/) || ($14 ~/^1\/1/ || $14 ~/^1\|1/) || ($15 ~/^1\/1/ || $15 ~/^1\|1/) || ($16 ~/^1\/1/ || $16 ~/^1\|1/) || ($17 ~/^1\/1/ || $17 ~/^1\|1/)  || ($18 ~/^1\/1/ || $18 ~/^1\|1/)  || ($19 ~/^1\/1/ || $19 ~/^1\|1/)  || ($20 ~/^1\/1/ || $20 ~/^1\|1/) || ($21 ~/^1\/1/ || $21 ~/^1\|1/) ){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t0|1"}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t1|0"}}' > BZ_NY_FixDiff.F1het.mm85.PHASE.vcf
Asks: if the BZ columns (columns 10-21 have 1/1 genotypes, meaning homozygous for the alt allele, then print the first 9 columns, and then "0|1", if the BZ columns do not have 1/1 (so they will be 0/0) then print "1|0" as the tenth column.
