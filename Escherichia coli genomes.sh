#!/bin/bash

# Initialization steps:
pushd ./alignment # Assuming the fastq.gz files are here
cp trimmomatic=$HOME/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa .
trimmomatic=$HOME/software/Trimmomatic-0.39/trimmomatic-0.39.jar
bowtie2-build ../sequence.fasta ecoli
samtools faidx ../sequence.fasta

# Run trimmomatic and bowtie for a given replicate
process_each(){
  # Assign name
  name=$1

  # Trim for quality
  java -jar $trimmomatic PE ${name}_FW.fastq.gz ${name}_RV.fastq.gz trim_${name}_FW_paired.fq.gz trim_${name}_FW_unpaired.fq.gz trim_${name}_RV_paired.fq.gz trim_${name}_RV_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

  # Align
  bowtie2 -x ecoli -1 trim_${name}_FW_paired.fq.gz -2 trim_${name}_RV_paired.fq.gz -S ${name}.sam

  # Convert to bam and sort
  samtools view -S -b ${name}.sam | samtools sort > ${name}.bam
}

# Call SNPs for a given replicate
call_each() {
  name=$1
  bcftools call -mv ${name}.bcf -Oz -o ${name}_call.vcf.bgz
  bcftools index ${name}_call.vcf.bgz
}

# Map Each
process_each control_rep_1
process_each control_rep_2
process_each control_rep_3

process_each resistant_rep_1
process_each resistant_rep_2
process_each resistant_rep_3

# samtools mpielup on each resistant .bam file
samtools mpileup -g -f ../sequence.fasta resistant_rep_1.bam > resistant_rep_1.bcf
samtools mpileup -g -f ../sequence.fasta resistant_rep_2.bam > resistant_rep_2.bcf
samtools mpileup -g -f ../sequence.fasta resistant_rep_3.bam > resistant_rep_3.bcf

# Call SNPs for each resistant .bcf file 
call_each resistant_rep_1
call_each resistant_rep_2
call_each resistant_rep_3

# Find SNPs present in all resistant replicates
bcftools isec -n=3 -c all -p dir resistant_rep_1_call.vcf.bgz resistant_rep_2_call.vcf.bgz resistant_rep_3_call.vcf.bgz

# samtools mpielup On each control .bam file
samtools mpileup -g -f ../sequence.fasta control_rep_1.bam > control_rep_1.bcf
samtools mpileup -g -f ../sequence.fasta control_rep_2.bam > control_rep_2.bcf
samtools mpileup -g -f ../sequence.fasta control_rep_3.bam > control_rep_3.bcf

# Call SNPs for each for each control .bcf file
call_each control_rep_1
call_each control_rep_2
call_each control_rep_3

# Subtract all control SNPs from resistant SNPs
bcftools view dir/0000.vcf | bgzip -c > dir/0000.vcf.bgz
bcftools index dir/0000.vcf.bgz
bcftools isec -n~1000 -c all -p dir_subtract dir/0000.vcf.bgz control_rep_1_call.vcf.bgz control_rep_2_call.vcf.bgz control_rep_3_call.vcf.bgz

