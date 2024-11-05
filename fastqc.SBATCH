#!/bin/bash
#SBATCH --job-name=fastqc.SBATCH--output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=jht47@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#---- This script runs fastqc on fq files ---------------------------------------#

#---- load module command -------------------------------------------------------#
module load fastqc

#---- DEFINE VARIABLES WITH FILE LOCATIONS --------------------------------------#
output_folder=/home/jht47/RNASeq_Project/fastqc_output
input=/home/jht47/RNASeq_Project/trimmomatic_output/WTB2_2.trPE.fq.gz

#---- RUN fastqc ----------------------------------------------------------------#
fastqc -o $output_folder $input

#---- unload module command -----------------------------------------------------#
module unload fastqc

#---- FIN -----------------------------------------------------------------------#
