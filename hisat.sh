#!/bin/bash
fastqc `ls -1 samples | cut -f1 -d"_" | uniq`
multiqc .
hisat2-build Bos_taurus.fa Bos_taurus.idx
for sample in `ls -1 samples | cut -f1 -d"_" | uniq`
do
 R1=samples/"${sample}"_1.fastq
 R2=samples/"${sample}"_2.fastq
 if [ -f $R1 ] && [ -f $R2 ]; then
  hisat2 -x Bos_taurus.idx -1 "${R1}" -2 "${R2}" -S "${sample}".sam --rg-id ${sample} --rg "PL:ILLUMINA" --rg SM:${sample}
 #Converting sam to bam
 samtools view -bS ${sample}.sam > ${sample}.bam
 #sorting bam files
 samtools sort ${sample}.bam > Bam_files/${sample}_sorted.bam
 samtools index ${sample}_sorted.bam
 gatk MarkDuplicates --INPUT ${sample}_sorted.bam --OUTPUT ${sample}_sorted_dedup.bam --METRICS_FILE  ${sample}_dup_metrics.txt --REMOVE_DUPLICATES true  --CREATE_INDEX true
  else
  echo "Error: $R1 or $R2 not found!"
  exit
  fi

done

gatk CreateSequenceDictionary -R Bos_taurus.fa
samtools faidx Bos_taurus.fa
wget https://ftp.ensembl.org/pub/release-108/variation/vcf/bos_taurus/bos_taurus.vcf.gz
gzip -d bos_taurus.vcf.gz
bgzip bos_taurus.vcf
tabix -f -p vcf bos_taurus.vcf.gz



for sample in `ls -1 samples | cut -f1 -d"_" | uniq`
do
    # adding N cigar reads
 gatk SplitNCigarReads -R Bos_taurus.fa -I ${sample}_sorted_dedup.bam -O ${Bam_files}_sorted_dedupCIGAR.bam
    # BaseRecalibration 
 gatk BaseRecalibrator --input ${sample}_sorted_dedupCIGAR.bam --output ${sample}_recal_data.table --reference Bos_taurus.fa --known-sites bos_taurus.vcf
 gatk ApplyBQSR --bqsr-recal-file ${sample}_recal_data.table --input ${sample}_sorted_dedupCIGAR.bam --output ${sample}_sorted_dedupCIGAR_BQSR_recal.bam --reference Bos_taurus.fa
 gatk BaseRecalibrator --input ${sample}_sorted_dedupCIGAR_BQSR_recal.bam --output ${sample}_post_recal_data.table --reference Bos_taurus.fa --known-sites bos_taurus.vcf
 gatk HaplotypeCaller --reference Bos_taurus.fa --input ${sample}_sorted_dedupCIGAR_BQSR_recal.bam --output ${sample}.g.vcf.gz --ERC GVCF
 gatk SelectVariants -R Bos_taurus.fa -V ${sample}_raw_variants.vcf --select-type-to-include SNP -O ${sample}_raw_snps.vcf
 gatk SelectVariants -R Bos_taurus.fa -V ${sample}_raw_variants.vcf --select-type-to-include INDEL -O ${sample}_raw_indels.vcf
 gatk VariantFiltration -V ${sample}_raw_snps.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter"FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum <-12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum <-8.0" --filter-name "ReadPosRankSum-8" -O ${sample}_snps_filtered.vcf
 gatk VariantFiltration -V ${sample}_raw_indels.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum <-20.0" --filter-name "ReadPosRankSum-20" -O ${sample}_indels_filtered.vcf
done

