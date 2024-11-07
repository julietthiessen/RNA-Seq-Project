#!/bin/bash
#SBATCH --job-name=trimmomatic.SBATCH --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=jht47@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#---- This script runs trimmomatic on paired end fastq files -------------#

#---- set environment, load module command -------------------------------#
shopt -s expand_aliases
module load trimmomatic

#---- DEFINE VARIABLES WITH FILE LOCATIONS -------------------------------------------------------#
adapters=/home/jht47/RNASeq_Project/TruSeq3-PE.fa

input_R1=/home/jht47/RNASeq_Project/WTB2_1.fq
input_R2=/home/jht47/RNASeq_Project/WTB2_2.fq

output_R1_PE=/home/jht47/RNASeq_Project/trimmomatic_output/WTB2_1.trPE.fq.gz
output_R1_SE=/home/jht47/RNASeq_Project/trimmomatic_output/WTB2_1.trSE.fq.gz

output_R2_PE=/home/jht47/RNASeq_Project/trimmomatic_output/WTB2_2.trPE.fq.gz
output_R2_SE=/home/jht47/RNASeq_Project/trimmomatic_output/WTB2_2.trSE.fq.gz

#---- RUN trimmomatic ----------------------------------------------------------------------------#
trimmomatic PE \
$input_R1 \
$input_R2 \
$output_R1_PE $output_R1_SE \
$output_R2_PE $output_R2_SE \
ILLUMINACLIP:$adapters:2:30:10 \
HEADCROP:15 \
TRAILING:20 \
SLIDINGWINDOW:4:15 \
MINLEN:75

#---- unload module command ----------------------------------------------------------------------#
module unload trimmomatic

#---- FIN ----------------------------------------------------------------------------------------#
