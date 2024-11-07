#!/bin/bash
#SBATCH --job-name=bowtie.SBATCH --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=jht47@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#---- This script runs bowtie2/2.5.3 on paired end fastq files -------------#

#---- load module command --------------------------------------------------#
module load bowtie2/2.5.3

#---- DEFINE VARIABLES WITH FILE LOCATIONS ---------------------------------#
index_prefix=/home/jht47/RNASeq_Project/bowtie_alignment/calb_ref_index

input_R1=/home/jht47/RNASeq_Project/trimmomatic_output/WTB2_1.trPE.fq
input_R2=/home/jht47/RNASeq_Project/trimmomatic_output/WTB2_2.trPE.fq

output=/home/jht47/RNASeq_Project/bowtie_alignment/WTB2.sam

#---- RUN bowtie -----------------------------------------------------------#
bowtie2 -x $index_prefix \
-1 $input_R1 \
-2 $input_R2 \
-S $output 

#---- unload module command ------------------------------------------------#
module unload bowtie2/2.5.3

#---- FIN ------------------------------------------------------------------#
