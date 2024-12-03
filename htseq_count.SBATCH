#!/bin/bash
#SBATCH --job-name=htseq_count.SBATCH --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=jht47@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#---- This script uses htseq to count and assign genes to reads ------------#

#---------------------------------------------------------------------------#
source activate htseq

#---- DEFINE VARIABLES WITH FILE LOCATIONS ---------------------------------#
htseq=/home/jht47/.conda/envs/htseq/bin/htseq-count

alignment_files=/home/jht47/RNASeq_Project/bowtie_alignment/WTB2.srt.bam
gtf_file=/home/jht47/RNASeq_Project/GCF_000182965.3_ASM18296v3_genomic.gtf

output=/home/jht47/RNASeq_Project/htseq_count/WTB2_htseqCount

#---- RUN htseq-count ------------------------------------------------------#
${htseq} -f bam -r pos -s no -t exon -i gene_id $alignment_files $gtf_file > $output

#---- FIN ------------------------------------------------------------------#
