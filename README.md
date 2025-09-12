# rup

RNA-seq usability assessment pipeline

Pipeline to evaluate the quality of RNA sequencing reads for differential expression analysis using different metrics. Only a minimal input file set is needed:

 - genome sequence in FASTA format
 - gene model annotation in GTF format
 - (optional) rRNA annotation in GTF format
 - sequencing reads in FASTQ.gz format
 <!-- - replicate to sample information in TSV format --> 

![RNA-seq usability pipeline](https://github.com/oliverrupp/rup/blob/main/images/Fig1.png?raw=true)


# Installation

The pipeline is implemented in R and depends on several R packages:

 - getopt (1.20.4, [CRAN](https://cran.r-project.org/web/packages/getopt/index.html))
 - ggplot2 (3.5.2, [CRAN](https://cran.r-project.org/web/packages/ggplot2/index.html))
 - reshape2 (1.4.4, [CRAN](https://cran.r-project.org/web/packages/reshape2/index.html))
 - pheatmap (1.0.13, [CRAN](https://cran.r-project.org/web/packages/pheatmap/index.html))
 - fastqcr (0.1.3, [CRAN](https://cran.r-project.org/web/packages/fastqcr/index.html))
 - Rfastp (1.16.0, [DOI](10.18129/B9.bioc.Rfastp))
 - Rsubread (2.20.0, [DOI](10.18129/B9.bioc.Rsubread))
 - Rsamtools (2.22.0, [DOI](10.18129/B9.bioc.Rsamtools))
 
All packages can be installed with [Bioconductor](https://bioconductor.org/): 

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("getopt", "ggplot2", "reshape2", "pheatmap", "fastqcr", "Rfastp", "Rsubread", "Rsamtools"))

fastqcr::fastqc_install()
```

Alternatively, dependencies can be install with conda:

```bash
conda env create -f conda.yaml
```


# Usage

All input files should be located relative to the data folder:


```bash
datafolder/
├── reads/
│   ├── s1_r1_1.fq.gz
│   ├── s1_r1_2.fq.gz
│   ├── s1_r2_1.fq.gz
│   ├── s1_r2_2.fq.gz
│   ├── s2_r1_1.fq.gz
│   ├── s2_r1_2.fq.gz
│   ├── s2_r2_1.fq.gz
│   ├── s2_r2_2.fq.gz
└── reference/
    ├── annotation.gtf
    ├── genome.fa
    └── rRNA.gtf
```

 - `reference/genome.fa` [genome FASTA file]
 - `reference/annotation.gtf` [gene model GTF file]
 - `reference/rRNA.gtf` [rRNA genes GTF file]
<!-- - `reference/samples.tsv` [replicate to sample mapping] -->
 - `reads/PREFIX_1.fastq.gz` and `reads/PREFIX_2.fastq.gz`


```
Usage: rup.R [-[-datafolder|d] <character>] [-[-threads|t] <integer>] [-[-bamSortMemory|m] <integer>] [-[-minReadLength|r] <integer>] [-[-minFragLength|l] <integer>] [-[-maxFragLength|u] <integer>] [-[-orientation|o] <character>] [-[-stranded|s] <integer>] [-[-help|h]]

Options:
 -d [folder]   location of the data folder
 -t [number]   number of threads to use
 -r [number]   minimum read length after trimming
 -m [number]   maximum memory for BAM file sorting
 -l [number]   minimum fragment length
 -u [number]   maximum fragment length
 -o [fr|rf|ff] paired end read orientation
 -h            print this help message
```

For single-end data use the `rup_SE.R` script. Each file in the `reads` folder will be analysed separately as single-end file.

```
Usage: rup_SE.R [-[-datafolder|d] <character>] [-[-threads|t] <integer>] [-[-bamSortMemory|m] <integer>] [-[-minReadLength|r] <integer>] [-[-orientation|o] <character>] [-[-stranded|s] <integer>] [-[-help|h]]

Options:
 -d [folder]   location of the data folder
 -t [number]   number of threads to use
 -r [number]   minimum read length after trimming
 -m [number]   maximum memory for BAM file sorting
 -o [fr|rf|ff] paired end read orientation
 -s [0|1|2]    stranded sequencing (0 (unstranded), 1 (stranded) and 2 (reversely stranded))
 -h            print this help message

```

# Prepare Files


## reference/genome.fa

The genome sequence must be provided in a single, uncompressed FASTA file.

## reference/annotation.gtf

The annotation should be in GTF format.
If only a GFF(3) format is available, the file can be converted to GTF using [gffread](https://github.com/gpertea/gffread).

```bash
gffread reference/annotation.gff3 -T -o reference/annotation.gtf
```

## reference/rRNA.gtf

rRNA genes can be predicted using [barrnap](https://github.com/tseemann/barrnap) and the provided `barrnap2gtf.pl` script.

```bash
barrnap --kingdom euk --threads 1 reference/genome.fa > reference/barrnap.gff

perl barrnap2gtf.pl reference/barrnap.gff > reference/rRNA.gtf
```

<!--
## reference/samples.tsv

A tab-delimited file with the sample names in the first column and the corresponding replicate names in the second column. The replicate names should be the same as the prefixes of the sequencing read files:

```bash
sample1<TAB>sample1_replicate1
sample1<TAB>sample1_replicate2
sample2<TAB>sample2_replicate1
sample2<TAB>sample2_replicate2
```
-->

## reads/PREFIX_[12].fastq.gz

The sequencing read files should be located in the `reads` folder.


# Results

All results and intermediate files will be saved in the `results` subfolder:

```bash
results/
├── trimmed/                      # FastP trimmed reads
├── fastqc/                       # FastQC analysis of the raw reads
├── trimmed_fastqc/               # FastQC analysis of the trimmed reads
├── bam/                          # unsorted read alignments in BAM format
├── sorted_bam/                   # sorted and indexed BAM files
└── counts/
    ├── feature.counts.tsv        # raw read counts
    ├── TPM.normalized.tsv        # TPM normalized read counts
    ├── sequencing_stats.tsv      # sequencing quality numbers
    ├── mapping_stats.tsv         # mapping quality numbers
    ├── gene_capture_stats.tsv    # number of captured genes
    └── sample_correlation.tsv    # pair-wise sample correlations
```


Additionally, the pipeline produces a PDF file (`RNAseq_QC.pdf`) in the data folder with the following analyses.

## sequencing quality 

![Number of reads before and after trimming](https://github.com/oliverrupp/rup/blob/main/images/Fig2.png?raw=true)

Assessment of the number of reads sequenced and remaining reads after trimming and filtering. 

## mapping quality

![Read Mapping Number](https://github.com/oliverrupp/rup/blob/main/images/Fig3.png?raw=true)

Assessment of the number of reads that could be assigned uniquely to a single gene.

## number of captured genes

![Number of Reads per Gene](https://github.com/oliverrupp/rup/blob/main/images/Fig4.png?raw=true)

Overview of the number of reads assigned to each gene.

## replicate quality

![Sample Correlation](https://github.com/oliverrupp/rup/blob/main/images/Fig5.png?raw=true)

Clustered heatmap of the replicate correlations


## RNA degradation

After running the pipeline, the included script `degradation_analysis_RSeQC.pl` can be used to check RNA degradation.
The script runs the `geneBody_coverage.py` and `tin.py` scripts from the [RSeQC package](https://rseqc.sourceforge.net/).

The path to the scripts must be included in the `$PATH` variable.

```bash
perl degradation_analysis_RSeQC.pl -d <DATAFOLDER>
```

The results will be stored in the `results/rseqc` subfolder.
