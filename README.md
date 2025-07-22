# rup

RNAseq usability pipeline

Pipeline to evaluate the quality of RNA sequencing reads for differential expression analysis using different metrics. Only a minimal input file set is needed:

 - genome sequence in FASTA format
 - gene model annotation in GTF format
 - (optional) rRNA annotation in GTF format
 - sequencing reads in FASTQ.gz format
 - replicate to sample information in TSV format

# Installation

The pipeline is implemented in R and depends on several R packages:

 - getopt ([CRAN](https://cran.r-project.org/web/packages/getopt/index.html))
 - ggplot2 ([CRAN](https://cran.r-project.org/web/packages/ggplot2/index.html))
 - reshape2 ([CRAN](https://cran.r-project.org/web/packages/reshape2/index.html))
 - pheatmap ([CRAN](https://cran.r-project.org/web/packages/pheatmap/index.html))
 - fastqcr ([CRAN](https://cran.r-project.org/web/packages/fastqcr/index.html))
 - Rfastp ([DOI](10.18129/B9.bioc.Rfastp))
 - Rsubread ([DOI](10.18129/B9.bioc.Rsubread))
 
All packages can be installed with [Bioconductor](https://bioconductor.org/): 

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("getopt", "ggplot2", "reshape2", "pheatmap", "fastqcr", "Rfastp", "Rsubread"))

fastqcr::fastqc_install()
```

# Usage

All input files should be located relative to the data folder:

 - reference/genome.fa [genome FASTA file]
 - reference/annotation.gtf [gene model GTF file]
 - reference/rRNA.gtf [rRNA genes GTF file]
 - reference/samples.tsv [replicate to sample mapping]
 - reads/PREFIX\_1.fastq.gz and reads/PREFIX\_2.fastq.gz

```bash
Usage: Rscript rup.R [-[-datafolder|d] <character>] [-[-threads|t] <integer>] [-[-help|h]]

Options:
 -d [folder] location of the data folder
 -t [number] number of threads to use
 -h          print this help message
```

# Prepare Files

## reference/genome.fa

The genome sequence must be provided in a single, uncompressed FASTA file.

## reference/annotation.gtf

The annotation should be in GTF format.
If only a GFF(3) format is available, the file can be converted to GTF using gffread

```bash
gffread reference/annotation.gff3 -T -o reference/annotation.gff
```

## reference/rRNA.gtf

rRNA genes can be predicted using barrnap and the provided barrnap2gtf.pl script.

```bash
barrnap --kingdom euk --threads 1 reference/genome.fa > reference/barrnap.gff

perl barrnap2gtf.pl reference/barrnap.gff > reference/rRNA.gtf
```

## reference/samples.tsv

A tab-delimited file with the sample names in the first column and the corresponding replicate names in the second column. The replicate names should be the same as the prefixes of the sequencing read files:

```bash
sample1<TAB>sample1_replicate1
sample1<TAB>sample1_replicate2
sample2<TAB>sample2_replicate1
sample2<TAB>sample2_replicate2
```

## reads/PREFIX_[12].fastq.gz

The sequencing read files should be located in the `reads` folder.


# Results

The pipeline produces a PDF file (`RNAseq_QC.pdf`) in the data folder with the following analyses.

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
