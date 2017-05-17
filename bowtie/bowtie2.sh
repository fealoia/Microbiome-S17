#!/bin/bash
#$ -l mem=10G,time=36:: -S /bin/bash -N bowtie2 -j y -cwd
#$ -o bowtie2.out
#$ -e bowtie2.err
bowtie2 --no-mixed --no-discordant --no-unal --no-hd --no-sq -x  $1 -1../../ip2169/Microbiome/poolA_S17_L005_R2_001.fastq -2../../ip2169/Microbiome/poolA_S17_L005_R1_001.fastq -S ${1}.sam
 
