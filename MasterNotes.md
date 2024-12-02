# RNA-Seq-Project

<h3>Goals of the Analysis:</h3>
<p>In the following analysis, we aim to characterize the distinct transcriptomes of populations of the commensal yeast <i>Candida albicans</i> grown in the presence and absence of vitamin B, or thiamine, using genomic analysis. <i>C. albicans</i> is a common cause of urinary tract, genital, mucosal, and blood infections in humans. Understanding how the yeast responds epigenetically to environmental changes via transcriptome modification is critical to evaluating and understanding the specific mechanisms behind potential treatment methods. Here, we will use bioinformatic programs to characterize the specific changes in gene expression associated with a thiamine-rich (Thi+) or thiamine-negative (Thi-) environments through differential expression analysis of the transcriptomes of <i>C. albicans</i> populations. This analysis will allow us to elucidate how transcription or expression of genes associated with the <i>C. albicans</i> response to thiamine changes in the presence of thiamine versus without thiamine, and identity which genes are specifically associated with the <i>C. albicans</i> thiamine response.</p>

<h3>Experimental Design:</h3>
<p>We obtained forward and reverse read files for six next-generation sequenced <i>C. albicans</i> RNA samples (.FQ format) from the Rolfes Laboratory. Three of the forward / reverse pair biological replicate files were obtained from samples or populations grown in the presence of thiamine and three of the forward / reverse pair biological replicate files were obtained from samples grown in the absence of thiamine.</p>
<p><b>File Names / Treatments:</b></p>
<p>WTA1_1.fq.gz, WTA1_2.fq.gz – C. albicans cultured in Thi+ environment; biological replicate A; forward and reverse reads</p>
<p>WTB1_1.fq.gz, WTB1_2.fq.gz – C. albicans cultured in Thi+ environment; biological replicate B; forward and reverse reads</p>
<p>WTC1_1.fq.gz, WTC1_2.fq.gz – C. albicans cultured in Thi+ environment; biological replicate C; forward and reverse reads</p>
<p>WTA2_1.fq.gz, WTA2_2.fq.gz – C. albicans cultured in Thi- environment; biological replicate A; forward and reverse reads</p>
<p>WTB2_1.fq.gz, WTB2_2.fq.gz – C. albicans cultured in Thi- environment; biological replicate B; forward and reverse reads</p>
<p>WTC2_1.fq.gz, WTC2_2.fq.gz – C. albicans cultured in Thi- environment; biological replicate C; forward and reverse reads</p>
<p>The following analysis was performed on the forward / reverse pair files for the second biological replicate data obtained from a <i>C. albicans</i> population grown in the absence of thiamine (forward read: WTB2_1; reverse read: WTB2_2) We ran the raw .FQ files through FastQC in order to determine the quality of the reads and developed an optimal trimming strategy. After trimming and clearing the reads using Trimmomatic, we confirmed improved quality again with FastQC. We then created an index and aligned the reads using Bowtie2/2.5.3 and Samtools/1.20 software. [Cont. - UNFINISHED]</p>

<h3>Data Files:</h3>
<p><b><i>C. albicans</i> sequenced RNA .FQ files were provided by the Rolfes Laboratory.</b></p>
<p><b><i>C. albicans</i> reference genome and annotation files (.FNA, .GTF) were obtained from the public NCBI database.</b></p>
<p>Link to database: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000182965.3/</p2>
<p>Reference genome assembly: ASM18296v3</p>
<p>NCBI RefSeq assembly: GCF_000182965.3</p>
<p>Link to downloadable files: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/</p>
<p>Full .FNA file name: GCF_000182965.3_ASM18296v3_genomic.fna.gz</p>
<p>Full .GTF file name: GCF_000182965.3_ASM18296v3_genomic.gtf.gz</p>

<h3>Workflow:</h3>
<p> 1. Quality Analysis (FastQC)</p>
<p>The quality of sequence reads for the obtained .FQ files were analyzed using FastQC software program on the Google Cloud HPC Console (bash). A trimming strategy was developed based on the results of the FastQC analysis.</p>
<p>Link to FastQC script used <a href=/fastqc.SBATCH> here</a>.</p>

<p> 2. Trimming & Quality Control (Trimmomatic)</p>
<p>A trimming strategy was developed according to the results of the FastQC analysis performed on the raw forward / reverse read files. After trimming was performed using Trimmomatic software on the Google Cloud HPC Console (bash), trimmed files were again analyzed using FastQC to verify improved quality of reads. The FastQC script used was identical to that provided above. The same Trimmomatic strategy was used on both the forward and the reverse reads.</p> 
<p>Link to Trimmomatic script used <a href=/trimmomatic.SBATCH> here<a/>.</p>

<p> 3. Alignment of Reads (Bowtie2/2.5.3 & Samtools/1.20)</p>
<p>A bowtie index file was created using the following command for Bowtie2/2.5.3 software on the Google Cloud HPC Console (bash) and the reference genome for <i>C. albicans</i> obtained from the public NCBI database.</p>
<p><b>Command: </b>$bowtie2-build GCF_000182965.3_ASM18296v3_genomic.fna calb_ref_index</p>
<p>Bowtie2/2.5.3 was then used to align the trimmed RNA sequence forward / reverse reads to the reference genome. The resulting output .SAM file was converted to .BAM format, sorted, and indexed using Samtools/1.20 software on the Google Cloud HPC Console (bash) using the following string of commands, respectively.</p>
<p><b>Commands: </b>$samtools view -S -b WTB2.sam > WTB2.bam; $samtools sort sample.bam –o sample.srt.bam; $samtools index sample.srt.bam<p/>
<p>Link to Bowtie2/2.5.3 script <a href=/bowtie.SBATCH> here<a/>.</p>

<p> 4. Counting of Reads Per Gene Model (HTseq)</p>
<p> 5. Differential Expression Analysis (DEseq2)</p>
<p> 6. Gene Ontology Enrichment Analysis</p>
 
<h3>Results:</h3>
<p> 1. Quality Analysis (FastQC)</p>
<p>Total sequences of forward read: 20353936</p>
<p>Total sequences of reverse read: 20353936</p>
<p>Issues identified: Per base sequence content, minor adaptor contamination</p>

<p> 2. Trimming & Quality Control (Trimmomatic)</p>
<p>Total sequences of forward read: 19476343</p>
<p>Total sequences of reverse read: 19476343</p>
<p>Issues identified: N/A</p>
<p>Link to Trimmomatic summary results: https://docs.google.com/spreadsheets/d/1AOa-XaTzR_PKMIRQDmu8oDTmawXXnkIwEjKOQkNC7Vs/edit?gid=0#gid=0</p>

<p> 3. Alignment of Reads (Bowtie2/2.5.3 & Samtools/1.20)</p>
<p>Total reads: 19476343</p>
<p>Number of read pairs detected: 19476343</p>
<p>Reads which aligned concordantly exactly 1 time: 89.04%</p>
<p>Reads which aligned concordantly > 1 time: 6.13%</p>
<p>Overall alignment rate: 98.11%</p>
<p>Link to full alignment results: https://docs.google.com/spreadsheets/d/1fa-FXVMlCXOZkbHSx_mMg0OXLMy9BeBJg8uWrEMpKGo/edit?gid=0#gid=0</p>

<p> 4. Counting of Reads Per Gene Model (HTseq)</p>

<p> 5. Differential Expression Analysis (DEseq2)</p>

<p> 6. Gene Ontology Enrichment Analysis</p>

<h3>Interpretation of Workflow:</h3>

<h3>Interpretation of Results:</h3>
