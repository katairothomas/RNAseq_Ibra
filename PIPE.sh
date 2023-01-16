#!/bin/bash
 for sample in `samples`;
 do
     R1=samples/"${sample}"_1.fastq
     R2=samples/"${sample}"_2.fastq
     #Alignment using Hisat2 while adding readgroups
  hisat2 -x Bos_taurus_index -1 "${sample}"_1.fastq -2 "${sample}"_2.fastq -S "${sample}".sam --rg-id RG1 --rg "PL:ILLUMINA" --rg SM:Sample1
  #Viewing the readgroups
  samtools view -H "${sample}".sam | grep '^@RG'
 #index for the genome file
 samtools faidx Bos_taurus.fa
 #sorting
 samtools sort "${sample}".bam > "${sample}"_sorted.bam
 #indexing sorted bam
 samtools index "${sample}"_sorted.bam
 #Creating a dictionary
 gatk CreateSequenceDictionary --REFERENCE Bos_taurus.fa --OUTPUT Bos_taurus.dict
 #Creating index file for vcf
 gatk IndexFeatureFile -I bos_taurus.vcf
 #marking duplicates
 gatk MarkDuplicates --INPUT "${sample}"_sorted.bam --OUTPUT gatk/"${sample}"_sorted_dedup.bam --METRICS_FILE metrics/"${sample}"_dup_metrics.txt --REMOVE_DUPLICATES true  --CREATE_INDEX true
 #creating N split cigar reads
 gatk SplitNCigarReads -R Bos_taurus.fa -I gatk/"${sample}"_sorted_dedup.bam -O "${sample}"CIGAR.bam
 #Base calibration
 gatk BaseRecalibrator --input gatk/"${sample}"_sorted_dedup.bam --output gatk/"${sample}"_recal_data.table --reference Bos_taurus.fa --known-sites bos_taurus.vcf
 #
 gatk ApplyBQSR --bqsr-recal-file gatk/"${sample}"_recal_data.table --input gatk/"${sample}"_sorted_dedup.bam --output gatk/"${sample}"_sorted_dedup_BQSR_recal.bam --reference Bos_taurus.fa
 gatk BaseRecalibrator --input gatk/"${sample}"_sorted_dedup_BQSR_recal.bam --output gatk/"${sample}"_post_recal_data.table --reference Bos_taurus.fa --known-sites bos_taurus.vcf
 #Haplotype caller
 gatk HaplotypeCaller --reference Bos_taurus.fa --input gatk/"${sample}"_sorted_dedup_BQSR_recal.bam --output gatk/"${sample}".g.vcf.gz --ERC GVCF

done
