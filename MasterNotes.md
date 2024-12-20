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
<p>The following analysis was performed on the forward / reverse pair files for the second biological replicate data obtained from a <i>C. albicans</i> population grown in the absence of thiamine (forward read: WTB2_1; reverse read: WTB2_2) We ran the raw .FQ files through FastQC in order to determine the quality of the reads and developed an optimal trimming strategy. After trimming and clearing the reads using Trimmomatic, we confirmed improved quality again with FastQC. We then created an index and aligned the reads using Bowtie2/2.5.3 and Samtools/1.20 software. HTSeq-Count was run on the resulting sorted, indexed read file to count the number of transcripts or reads per gene model. The “WTB2” htseqCount output file was then recombined with the remaining five files and the following analysis was performed on all six files together. In order to identify the differentially expressed genes in the counted files, we used an R-script to run DESeq2, produce a PCA plot, produce a volcano plot, and identify differentially expressed genes. We analyzed the biological function of the identified genes both manually, using the NCBI and Candida Genome databases, and quantitatively, by performing a Gene Ontology Enrichment Analysis on the Panther database website.</p>

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
<p>Link to Trimmomatic summary results: https://docs.google.com/spreadsheets/d/1AOa-XaTzR_PKMIRQDmu8oDTmawXXnkIwEjKOQkNC7Vs/edit?gid=0#gid=0</p>
<p>New file names: WTB2_1.trPE.fq; WTB2_2.trPE.fq</p>

<p> 3. Alignment of Reads (Bowtie2/2.5.3 & Samtools/1.20)</p>
<p>A bowtie index file was created using the following command for Bowtie2/2.5.3 software on the Google Cloud HPC Console (bash) and the reference genome for <i>C. albicans</i> obtained from the public NCBI database.</p>
<p><b>Command: </b>$bowtie2-build GCF_000182965.3_ASM18296v3_genomic.fna calb_ref_index</p>
<p>Bowtie2/2.5.3 was then used to align the trimmed RNA sequence forward / reverse reads to the reference genome. The resulting output .SAM file was converted to .BAM format, sorted, and indexed using Samtools/1.20 software on the Google Cloud HPC Console (bash) using the following string of commands, respectively.</p>
<p><b>Commands: </b>$samtools view -S -b WTB2.sam > WTB2.bam; $samtools sort sample.bam –o sample.srt.bam; $samtools index sample.srt.bam<p/>
<p>Link to Bowtie2/2.5.3 script <a href=/bowtie.SBATCH> here<a/>.</p>
<p>Link to full alignment results: https://docs.google.com/spreadsheets/d/1fa-FXVMlCXOZkbHSx_mMg0OXLMy9BeBJg8uWrEMpKGo/edit?gid=0#gid=0</p>
<p>New file names: WTB2.sam; WTB2.bam; WTB2.srt.bam; WTB2.srt.bam.bai</p>

<p> 4. Counting of Reads Per Gene Model (htseqCount)</p>
<p>Conda software was used to install the compartmentalized suite of softwares to structure the environment for the HTSeq python program using the bioconda package “htseq.” The htseq-count command was run on the Google Cloud HPC Console (bash) using the following script to parse through the WTB2.str.bam file, match the reads with gene locations in the annotated C. albicans .gtf file, and count the number of transcripts or reads per gene model. The htseqCount outputs files for all six original read files were combined in a singular folder and used together from this point on in the workflow.</p>
<p>Link to HTSeq script <a href=/htseq_count.SBATCH> here<a/>.</p>
<p>New files: NTH_WTA2_htseqCount; NTH_WTB2_htseqCount; NTH_WTC2_htseqCount; TH_WTA1_htseqCount; TH_WTB1_htseqCount; TH_WTC1_htseqCount</p>

<p> 5. Differential Expression Analysis (DEseq2)</p>
<p>The differential expression analysis program (DESeq2) was run using an R-script on the R Studio platform. A DESeq data set was generated using the input HTSeq Count data and associated files. The data was then filtered and provided a reference condition (the control data, from yeast grown in the presence of thiamine). Quality control was performed on the data, and a PCA visualization plot was generated for the 6 samples (see results). The Differentially Expressed Genes analysis was then performed using a negative binomial distribution to model the RNA-Seq counts. Output genes in the table were ranked by smallest p-value, and written into a table (see results). The differentially expressed genes were mapped on a volcano plot (see results). Finally, gene names were added to the differentially expressed genes according to information pulled from the NCBI database and a summary table of results was generated (see results).</p>
<p>Link to DESeq2 script <a href=/calb_DESeq_script_FINAL.R> here<a/>.</p>
<p>Finally, the resulting summary table was manually supplemented with the NCBI Gene ID, GTF annotation file gene name and a description of the biological function of the gene. To obtain the gene ID numbers and annotation file gene names, a file “signif_geneIDs” was created on the Google Cloud HPC Console (bash) with the output Locus Tags provided for the differentially expressed genes. The following bash command and GCF_000182965.3_ASM18296v3_genomic.gtf.gz annotation file was then used to generate a new file “signif_gene_annot_info” with the desired information. Finally, information from the NCBI database and online Candida Genome database was used to supplement the table with the biological functions of the identified genes.</p>
<p><b>Command: </b>$grep -wFf signif_geneIDs GCF_000182965.3_ASM18296v3_genomic.gtf|grep
"protein_coding"|cut -f9|cut -d ";" -f1,3,5 > signif_gene_annot_info</p>
<p>NCBI database: https://www.ncbi.nlm.nih.gov/search/</p>
<p>Candida Genome database: http://www.candidagenome.org/</p>
 
<p> 6. Gene Ontology Enrichment Analysis</p>
<p>The Panther classification system website / database was used to perform a Gene Ontology Enrichment Analysis on the 13 differentially expressed genes identified using DESeq2 (above). NCBI GeneID #s for each of the genes were input, the organism <i>Candida albicans</i> was selected and a Statistical Overrepresentation test was run to identify overrepresented or enriched (biological process) GO terms for the set of differentially expressed genes. “GO biological process complete” was selected to identify terms relevant to the biological function of the genes. The <i>C. albicans</i> complete genome set was selected as a reference and a statistical Fisher test calculating a false discovery rate was performed. The resulting enriched terms were downloaded into an Excel file for analysis (see results).</p>
<p>Panther website: http://www.pantherdb.org/</p>
 
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

<p> 4. Counting of Reads Per Gene Model (htseqCount)</p>
<p>No feature: 369641</p>
<p>Ambiguous: 286389</p>
<p>Too low quality: 1083591</p>
<p>Not aligned: 165060</p>
<p>Alignment not unique: 0</p>

<p> 5. Differential Expression Analysis (DEseq2)</p>
<p>Link to PCA plot <a href=/TH-vTH+_pcaplot.png> here<a/>.</p>
<p>Link to complete gene table <a href=/calb_TH-_v_TH+_allgenes.csv> here<a/>.</p>
<p>Link to volcano plot of differentially expressed genes <a href=/R_volcano_plot_correct.pdf> here<a/>.</p>
<p>Link to table of differentially expressed genes <a href=/signif_TH-vTH+.csv> here<a/>.</p>
<p>Link to supplemented table of differentially expressed genes with function <a href=/Summary_Table_(signif_TH-vTH+).xls> here<a/>.</p>

<p> 6. Gene Ontology Enrichment Analysis</p>
<p>Link to GO enrichment analysis Panther summary table <a href=/GO_Enrichment_Analysis_Results.xls> here<a/>.</p>

<h3>Interpretation of Results:</h3>
<p> 5. Differential Expression Analysis (DEseq2)</p>
<p>We identified several patterns from the identification of the biological function of significantly differentially expressed genes (Workflow / Results 5). Of the thirteen significantly differentially expressed genes identified in <i>Candida albicans</i> grown in thiamine-negative environments, functional information was available for nine of the genes. Taken together, patterns were identified across the functional descriptions for the nine genes (see results). First, seven of the nine genes are known to be involved in Vitamin B biosynthetic processes and included terms such as thiamine synthesis, pyridoxine biosynthesis and metabolism, Vitamin B synthesis, thiamine biosynthesis, hydroxymethylpyrimidine phosphate thiamine precursor, thiamine-phosphate pyrophosphorylase, pyridoxal phosphate biosynthesis, and hydroxyethylthiazole kinase activity. Each of these processes are related to the biosynthesis, pathway, and metabolism of various forms of Vitamin B. When grown in a thiamine-negative environment (lacking Vitamin B), strains of <i>C. albicans</i> upregulate the production of proteins involved in the synthesis and metabolism of Vitamin B forms internally which do not need to be expressed when the yeast is grown in a thiamine-rich environment. Second, we identified one significantly differentially expressed gene involved in Vitamin B transmembrane transport into the cell (FCY24). This discovery indicates that the <i>C. albicans</i> grown in thiamine-negative environments may be upregulating Vitamin B transporters in an attempt to account for an environment which may have less available thiamine and thus required additional transporters. Third and finally, seven of the nine genes identified were “spider biofilm induced.” <i>C. albicans</i> yeast grown in the nutrient-deprived thiamine-negative environment also responded to the environmental change by upregulating spider biofilm as a defensive mechanism and to induce the above expression changes.</p>

<p> 6. Gene Ontology Enrichment Analysis</p>
<p>The Gene Ontology Enrichment Analysis performed for the identified differentially expressed genes further validated the above conclusions. We performed the GO Enrichment Analysis using the Panther website for the Biological Processes category to identify GO terms associated with the functions of the proteins produced by the differentially expressed genes identified for <i>C. albicans</i> grown in a thiamine-negative environment. Within the thirteen genes, there were more instances of these associated terms than would be expected in a random sample of thirteen genes. 37 total Biological Processes GO terms were identified as “enriched” in the thirteen differentially expressed genes (see results). The terms were grouped into four subcategories: thiamine biosynthetic process, pyridoxine biosynthetic process, pyridoxal phosphate biosynthetic process, and ergosterol biosynthetic process. Thiamine, pyridoxine, and pyridoxal phosphate are associated with internal production of Vitamins B1 and B6; yeast grown in a thiamine-negative environment must produce Vitamin B themselves. Ergosterol is a critical protein in yeast cell membranes which is involved in the transport of proteins from the environment into the cell; <i>C. albicans</i> grown in thiamine-negative environments increase production of transmembrane proteins in order to upregulate import of Vitamin B in the extracellular environment.</p>
