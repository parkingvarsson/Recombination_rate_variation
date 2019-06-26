#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00

module load bioinfo-tools
module load samtools/1.6

##### Script for splitting the fasta file into chromosome fasta files #####
samtools faidx F1.chr.fasta

samtools faidx F1.chr.fasta chr1 > chr1.fa
samtools faidx F1.chr.fasta chr2 > chr2.fa
samtools faidx F1.chr.fasta chr3 > chr3.fa
samtools faidx F1.chr.fasta chr4 > chr4.fa
samtools faidx F1.chr.fasta chr5 > chr5.fa
samtools faidx F1.chr.fasta chr6 > chr6.fa
samtools faidx F1.chr.fasta chr7 > chr7.fa
samtools faidx F1.chr.fasta chr8 > chr8.fa
samtools faidx F1.chr.fasta chr9 > chr9.fa
samtools faidx F1.chr.fasta chr10 > chr10.fa
samtools faidx F1.chr.fasta chr11 > chr11.fa
samtools faidx F1.chr.fasta chr12 > chr12.fa
samtools faidx F1.chr.fasta chr13 > chr13.fa
samtools faidx F1.chr.fasta chr14 > chr14.fa
samtools faidx F1.chr.fasta chr15 > chr15.fa
samtools faidx F1.chr.fasta chr16 > chr16.fa
samtools faidx F1.chr.fasta chr17 > chr17.fa
samtools faidx F1.chr.fasta chr18 > chr18.fa
samtools faidx F1.chr.fasta chr19 > chr19.fa
