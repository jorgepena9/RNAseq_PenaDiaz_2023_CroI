# RNA-seq analysis - Alignment
#### Jorge Peña-Díaz - May 14th 2023

-----------------------------------------

Analysis used to compare the gene expression profile of the ∆*croI* and wild-type (WT) *C. rodentium* DB1000 strains as published in:

> Peña‑Díaz, J., Woodward, S.E., Creus‑Cuadros, A., Serapio‑Palacios, A., Deng, W., and Finlay, B.B. (2023). Quorum Sensing Modulates
Bacterial Virulence and Colonization Dynamics During an Enteric Infection. bioRxiv, [10.1101/2023.03.14.532656](https://www.biorxiv.org/content/10.1101/2023.03.14.532656v1).


-----------------------------------------

### 1. Environment

Define a global working directory and other relevant environment variables used throughout this analysis. For easier analysis, you can add these environmental variables to your `.bashrc` file.

```
export RNA_HOME=~/scratch/rnaseq_citrobacter
export RNA_DATA_DIR=$RNA_HOME/data
export RNA_REFS_DIR=$RNA_HOME/refs
export RNA_REFS_DIR_ICC168=$RNA_REFS_DIR/ncbi_dataset/data/GCF_000027085.1
export RNA_REFS_DIR_DBS100=$RNA_REFS_DIR/ncbi_dataset/data/GCF_000835925.1
export RNA_REFS_DIR_NBRC105723=$RNA_REFS_DIR/ncbi_dataset/data/GCF_000759815.1
export RNA_DATA_TRIM_DIR=$RNA_DATA_DIR/trimmed
export RNA_ALIGN_DIR=$RNA_HOME/alignments/STAR
```


### 2. Create working directory and raw data

This workflow starts with the raw paired-end data in demultiplexed FASTQ format assumed to be located within a folder in `$RNA_HOME/raw_data`. The RNA-seq data generated for this study is available on the NCBI Sequence Read Archive (SRA) under the BioProject accession number PRJNA935367. 

```
# Creating a working directory
mkdir -p $RNA_HOME
cd $RNA_HOME

# Create data folder
mkdir -p $RNA_DATA_DIR
cd $RNA_DATA_DIR
```

### 3. Pre-alignment QC

Inspect read quality before alignment

```
# Go to the data folder
cd $RNA_DATA_DIR

# Running FastQC
fastqc -t 8 *.fastq.gz

# Run MultiQC
multiqc .

```

### 4. Reference genome 

Download the *C. rodentium* ICC168 reference genome by utilizing [NCBI datasests](https://www.ncbi.nlm.nih.gov/datasets/).

```
# Create the reference genome folder
mkdir $RNA_REFS_DIR
cd $RNA_REFS_DIR

# Download Citrobacter genomes from NCBI
datasets download genome taxon 67825 \
--assembly-source refseq \
--annotated \
--assembly-level complete_genome,contig \
--include-gtf \
--filename citrobacter-rodentium.zip

# Unzip the data
unzip citrobacter-rodentium.zip

```

### 5. Alignment

Perform read alignment by using STAR. 

**Note:** As we are working with bacteria and their transcripts consist of single exons without splices, we replaced all features in the third column of the `.gtf` file with "exon", and then used the command `--sjdbGTFfeatureExon exon` during genome indexing.

```
# Move to required folder
cd $RNA_HOME

# Genome Index
STAR --runMode genomeGenerate \
--genomeDir $RNA_REFS_DIR_ICC168 \
--genomeFastaFiles $RNA_REFS_DIR_ICC168/GCF_000027085.1_ASM2708v1_genomic.fna \
--sjdbGTFfile $RNA_REFS_DIR_ICC168/genomic_exon.gtf \
--runThreadN 16 \
--alignIntronMax 1 \
--genomeSAindexNbases 8 \
--sjdbGTFfeatureExon exon

# Read alignment
for i in $RNA_DATA_DIR/*R1.fastq.gz ; do name=$(basename ${i} _R1.fastq.gz); name2=$(basename ${i} R1.fastq.gz | cut -b 8- | rev | cut -c 6- | rev);
STAR --genomeDir $RNA_REFS_DIR_ICC168/ \
--quantMode GeneCounts \
--runThreadN 16 \
--limitBAMsortRAM 20000000000 \
--readFilesIn $RNA_DATA_DIR/${name}_R1.fastq.gz,$RNA_DATA_DIR/January_2022/${name2}_R1.fastq.gz $RNA_DATA_DIR/${name}_R2.fastq.gz,$RNA_DATA_DIR/January_2022/${name2}_R2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix $RNA_ALIGN_DIR/ICC168/${name2} \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes Standard \
--twopassMode Basic ; done

```

### 5. Alignment QC

Merge, and index `.bam` files for alignment visualization

```
# Move to required folder 
cd $RNA_ALIGN_DIR/ICC168

# Create required folder
mkdir -p $RNA_HOME/results/merged_bam/
  
# Merge .bam files
for i in *_1Aligned.sortedByCoord.out.bam; do name=$(basename ${i} _1Aligned.sortedByCoord.out.bam);
java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=$RNA_HOME/results/merged_bam/${name}.bam INPUT=${name}_1Aligned.sortedByCoord.out.bam INPUT=${name}_3Aligned.sortedByCoord.out.bam INPUT=$RNA_ALIGN_DIR/ICC168/${name}_5Aligned.sortedByCoord.out.bam; done


# Index all .bam files
cd $RNA_HOME/results/merged_bam/
find *.bam -exec echo samtools index {} \; | sh
cd $RNA_ALIGN_DIR/ICC168
find *.bam -exec echo samtools index {} \; | sh

```

In order to generate RNA-seq specific quality metrics, we must first create the required files.


```
# For Picard

# Move to required folder 
cd $RNA_REFS_DIR_ICC168/
  
# Create a .dict file for our reference
java -jar $PICARD CreateSequenceDictionary R=GCF_000027085.1_ASM2708v1_genomic.fna O=GCF_000027085.1_ASM2708v1_genomic.dict

# Create a bed file of the location of ribosomal sequences in our reference (first extract from the gtf then convert to bed)
grep --color=none -i rrna genomic.gtf > ref_ribosome.gtf
gff2bed < ref_ribosome.gtf > ref_ribosome.bed

# Create interval list file for the location of ribosomal sequences in our reference
java -jar $PICARD BedToIntervalList I=ref_ribosome.bed O=ref_ribosome.interval_list SD=GCF_000027085.1_ASM2708v1_genomic.dict

# Create a genePred file for our reference transcriptome
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons genomic.gtf genomic.ref_flat.txt

# reformat this genePred file
cat genomic.ref_flat.txt | awk '{print $12"\t"$0}' | cut -d$'\t' -f1-11 > tmp.txt
mv tmp.txt genomic.ref_flat.txt


# For RSeQC

# Move to required folder 
cd $RNA_REFS_DIR_ICC168/

# Convert Gtf to genePred
gtfToGenePred -ignoreGroupsWithoutExons genomic.gtf genomic.genePred

# Convert genPred to bed12
genePredToBed genomic.genePred genomic.bed12

```








Use FastQC, Flagstat, Picard, RSeQC and MultiQC to evaluate the quality of the alignments of the `.bam` files.

```
# FastQC
fastqc -t 16 *.bam


# Flagstat
mkdir flagstat
find .bam -exec echo samtools flagstat {} \> flagstat/{}.flagstat \; | sh


# Picard - CollectRnaSeqMetrics
mkdir picard
find *.bam -exec echo java -jar $PICARD CollectRnaSeqMetrics I={} O=picard/{}.RNA_Metrics REF_FLAT=$RNA_REFS_DIR_ICC168/genomic.ref_flat.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=$RNA_REFS_DIR_ICC168/ref_ribosome.interval_list \; | sh


# RSeQC
mkdir rseqc
find *.bam -exec echo bam_stat.py -i {} \> {}.bam_stat.txt \; | sh
find *.bam -exec echo infer_experiment.py -i {} -r $RNA_REFS_DIR_ICC168/genomic.bed12 \> rseqc/{}.strandness.txt \; | sh
find *.bam -exec echo read_distribution.py  -i {} -r $RNA_REFS_DIR_ICC168/genomic.bed12 \> rseqc/{}.read_dist.txt \; | sh

# (Optional: These are very computationally intensive)
# for i in *_1Aligned.sortedByCoord.out.bam; do name=$(basename ${i} _1Aligned.sortedByCoord.out.bam); 
# geneBody_coverage.py -i ${name}_1Aligned.sortedByCoord.out.bam,${name}_3Aligned.sortedByCoord.out.bam,${name}_5Aligned.sortedByCoord.out.bam -r $RNA_REFS_DIR_ICC168/genomic.bed12 -o rseqc/${name}; done


# MultiQC
multiqc .

```

Gene counts were then imported to R for further processing.

-----------------------------------------

### 6. Software

* fastqc/0.11.9

* python/3.10.2

* multiQC/1.12

* datasets/13.17.0

* star/2.7.9a

* picard/2.23.3

* samtools/1.12

* bedops/2.4.39

* kentutils/401

* bedtools/2.29.2

* RSeQC/4.0.0

-----------------------------------------

### 7. Acknowledgements
This analysis assumes the use of a Linux computer with an ‘x86_64’ architecture. This research was partly enabled by software provided by the  [Digital Research Alliance of Canada](https://alliancecan.ca/en).

Sections of this pipeline were adapted from:

* Griffith M, Walker JR, Spies NC, Ainscough BJ, Griffith OL (2015) Informatics for RNA Sequencing: A Web Resource for Analysis on the Cloud. PLOS Computational Biology 11(8): e1004393. [https://doi.org/10.1371/journal.pcbi.1004393](https://doi.org/10.1371/journal.pcbi.1004393)

* Mary E. Piper, Meeta Mistry, Jihe Liu, William J. Gammerdinger, & Radhika S. Khetani. (2022, January 10). hbctraining/Intro-to-rnaseq-hpc-salmon-flipped: Introduction to RNA-seq using Salmon Lessons from HCBC (first release). Zenodo. [https://doi.org/10.5281/zenodo.5833880](https://doi.org/10.5281/zenodo.5833880)

