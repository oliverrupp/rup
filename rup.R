library(getopt)

spec = matrix(c(
  'datafolder' , 'd', 1, "character",
  'threads', 't', 1, "integer",
  'bamSortMemory', 'm', 1, 'integer',
  'minReadLength', 'r', 1, 'integer',
  'minFragLength', 'l', 1, 'integer',
  'maxFragLength', 'u', 1, 'integer',
  'orientation', 'o', 1, 'character',
  'stranded', 's', 1, 'integer',
  'help', 'h', 0, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

source_folder <- "."
n_threads     <- 1
bamSortMemory <- "1024"
minReadLength <- 25
minFragLength <- 0
maxFragLength <- 300
orientation   <- "fr"
stranded      <- 0

if (!is.null(opt$threads)) n_threads <- opt$threads
if (!is.null(opt$datafolder)) setwd(opt$datafolder)
if (!is.null(opt$orientation)) orientation <- opt$orientation
if (!is.null(opt$stranded)) stranded <- opt$stranded
if (!is.null(opt$maxFragLength)) maxFragLength <- opt$maxFragLength
if (!is.null(opt$minFragLength)) minFragLength <- opt$minFragLength
if (!is.null(opt$minReadLength)) minReadLength <- opt$minReadLength
if (!is.null(opt$bamSortMemory)) bamSortMemory <- opt$bamSortMemory

reference_folder      <- file.path(source_folder, "reference")
read_file_folder      <- file.path(source_folder, "reads")
results_folder        <- file.path(source_folder, "results")

genome_fasta_file     <- file.path(reference_folder, "genome.fa")
annotation_file       <- file.path(reference_folder, "annotation.gtf")
rrna_file             <- file.path(reference_folder, "rRNA.gtf")
### samples_file          <- file.path(reference_folder, "samples.txt")

fastqc_folder         <- file.path(results_folder,   "fastqc")
trimmed_fastqc_folder <- file.path(results_folder,   "trimmed_fastqc")
trimmed_read_folder   <- file.path(results_folder,   "trimmed")
bam_folder            <- file.path(results_folder,   "bam")
sorted_bam_folder     <- file.path(results_folder,   "sorted_bam")
counts_folder         <- file.path(results_folder,   "counts")

sample_prefixes       <- gsub("_1.f(ast)?q.gz", "", (list.files(read_file_folder, pattern = "*_1.f(ast)?q.gz")))

error <- FALSE
help  <- !is.null(opt$help)

if(!help && !file.exists(genome_fasta_file)) {
    message(paste0("ERROR: genome FASTA file is file missing [",genome_fasta_file,"]!"))
    error <- TRUE
}

if(!help && !file.exists(annotation_file)) {
    message(paste0("ERROR: annotation GTF file is file missing [",annotation_file,"]!"))
    error <- TRUE
}

if(!help && length(sample_prefixes) == 0) {
    message(paste0("ERROR: no FASTQ file found in ",read_file_folder,"!"))
    error <- TRUE
}

if(error || help) {
    message()
    cat(getopt(spec, usage=TRUE))
    message()
    message("Options:")
    message(" -d [folder]   location of the data folder")
    message(" -t [number]   number of threads to use")
    message(" -r [number]   minimum read length after trimming")
    message(" -m [number]   maximum memory for BAM file sorting")
    message(" -l [number]   minimum fragment length")
    message(" -u [number]   maximum fragment length")
    message(" -o [fr|rf|ff] paired end read orientation")
    message(" -s [0|1|2]    stranded sequencing (0 (unstranded), 1 (stranded) and 2 (reversely stranded))")
    message(" -h            print this help message")
    message()
    q(status=1)
}

for(folder in c(results_folder, fastqc_folder, trimmed_fastqc_folder, trimmed_read_folder, bam_folder, sorted_bam_folder, counts_folder)) {
  if(!dir.exists(folder)) {
    dir.create(folder)
  }
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(Rsamtools))

library(pheatmap)
library(fastqcr)
library(Rfastp)
library(Rsubread)

theme_set(theme_bw())



pdf("RNAseq_QC.pdf", w=18, h=12)


#########################################################################################################
#
# read based statistics
#
#########################################################################################################

#########################################################################################################
#### fastqc raw data

message("running fastqc ...")
fastqc_raw.start <- Sys.time()
fastqc(fq.dir=read_file_folder, qc.dir=fastqc_folder, threads=n_threads)

qc <- qc_aggregate(fastqc_folder, progress=F)
qc$tot.seq <- as.numeric(qc$tot.seq)
qc$sample = gsub(".f(ast)?q.gz", "" , qc$sample)

fastqc_raw.end <- Sys.time()
cat(sprintf("FastQC Raw: %0.2f minutes\n", difftime(fastqc_raw.end, fastqc_raw.start, units="mins")), file="benchmark.log")

#########################################################################################################
#### trimming

trimming.start <- Sys.time()
for(prefix in sample_prefixes) {
  outputPrefix <- file.path(trimmed_read_folder, prefix)

  if(!file.exists(paste(outputPrefix,"_R1.fastq.gz", sep=""))) {

    sample_trim.start <- Sys.time()
    read1 = file.path(read_file_folder, paste(prefix, "_1.fastq.gz", sep=""))
    read2 = file.path(read_file_folder, paste(prefix, "_2.fastq.gz", sep=""))
    if(!file.exists(read1)) { read1 = file.path(read_file_folder, paste(prefix, "_1.fq.gz", sep="")) }
    if(!file.exists(read2)) { read2 = file.path(read_file_folder, paste(prefix, "_2.fq.gz", sep="")) }
    
    fastp_stats <- rfastp(read1 = read1,
                          read2 = read2,
			  minReadLength = minReadLength,
                          outputFastq = outputPrefix, 
			  thread = n_threads)
    sample_trim.end <- Sys.time()
    cat(sprintf("Trimming (%s): %0.2f minutes\n", prefix, difftime(sample_trim.end, sample_trim.start, units="mins")), file="benchmark.log", append=T)
  }
}
trimming.end <- Sys.time()
cat(sprintf("Trimming: %0.2f minutes\n", difftime(trimming.end, trimming.start, units="mins")), file="benchmark.log", append=T)


#########################################################################################################
#### fastqc trimmed data

fastqc_trimmed.start <- Sys.time()
fastqc(fq.dir=trimmed_read_folder, qc.dir=trimmed_fastqc_folder, threads=n_threads)

qc_trimmed <- qc_aggregate(trimmed_fastqc_folder, progress=F)
qc_trimmed$tot.seq <- as.numeric(qc_trimmed$tot.seq)
qc_trimmed$sample = gsub("_R([12])$", "_\\1" , qc_trimmed$sample)

qc$Trimming = "raw"
qc_trimmed$Trimming = "trimmed"

qc_all = rbind(qc, qc_trimmed)

read_number_plot <- ggplot(qc_all, aes(x=sample, y=tot.seq, fill=Trimming)) + 
  geom_bar(stat="identity", position = "dodge") + 
  xlab("Samples") + ylab("Number of Reads") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  theme(text = element_text(size = 18)) + 
  ggtitle("Number of reads before and after trimming")
print(read_number_plot)

write.table(qc_all, "results/counts/sequencing_stats.tsv", quote=F, sep="\t", row.names=F)

fastqc_trimmed.end <- Sys.time()
cat(sprintf("FastQC Trimmed: %0.2f minutes\n", difftime(fastqc_trimmed.end, fastqc_trimmed.start, units="mins")), file="benchmark.log", append=T)

#########################################################################################################
#
# mapping based statistics
#
#########################################################################################################

#########################################################################################################
### read alignment

index.start <- Sys.time()
if(!file.exists(file.path(reference_folder,"subread.index.reads"))) {
  buildindex(file.path(reference_folder,"subread.index"), genome_fasta_file)
}
index.end <- Sys.time()
cat(sprintf("Genome index: %0.2f minutes\n", difftime(index.end, index.start, units="mins")), file="benchmark.log", append=T)

mapping.start <- Sys.time()
for(prefix in sample_prefixes) {
  message(prefix)

  output_bam = file.path(bam_folder, paste(prefix, ".bam", sep=""))
  output_sorted_bam = file.path(sorted_bam_folder, prefix)

  if(!file.exists(output_bam)) {
    sample_map.start <- Sys.time()
    mapping_stats <- align(index=file.path(reference_folder,"subread.index"),
                           readfile1=file.path(trimmed_read_folder, paste(prefix, "_R1.fastq.gz", sep="")),
                           readfile2=file.path(trimmed_read_folder, paste(prefix, "_R2.fastq.gz", sep="")),
                           output_file=output_bam, 
			   type=0,
                           minFragLength=minFragLength, 
			   maxFragLength=maxFragLength, 
			   PE_orientation=orientation,
                           nthreads=n_threads,
                           useAnnotation=TRUE,
                           annot.ext=annotation_file,
                           isGTF=TRUE,
                           nBestLocations=2)
    sample_map.end <- Sys.time()
    cat(sprintf("Mapping (%s): %0.2f minutes\n", prefix, difftime(sample_map.end, sample_map.start, units="mins")), file="benchmark.log", append=T)
  }

  if(!file.exists(paste0(output_sorted_bam, ".bam"))) {
	sample_sort.start <- Sys.time()
	sortBam(output_bam, output_sorted_bam, maxMemory=bamSortMemory, nThreads=n_threads)
	indexBam(paste0(output_sorted_bam, ".bam"))
	sample_sort.end <- Sys.time()
    	cat(sprintf("Sorting (%s): %0.2f minutes\n", prefix, difftime(sample_sort.end, sample_sort.start, units="mins")), file="benchmark.log", append=T)
  }
}
mapping.end <- Sys.time()
cat(sprintf("Mapping: %0.2f minutes\n", difftime(mapping.end, mapping.start, units="mins")), file="benchmark.log", append=T)




#########################################################################################################
### read counting

bam_file <- list.files(bam_folder, pattern = "*.bam$")

featurecounts.start <- Sys.time()
gene_feature_counts = featureCounts(file.path(bam_folder, bam_file),
                                    annot.ext              = annotation_file,
                                    isGTFAnnotationFile    = TRUE,
                                    countMultiMappingReads = FALSE,
                                    strandSpecific         = stranded,
                                    isPairedEnd            = TRUE,
				    requireBothEndsMapped  = TRUE,
				    checkFragLength        = TRUE,
				    minFragLength          = minFragLength,
				    maxFragLength          = maxFragLength,
                                    nthreads               = n_threads)

featurecounts.end <- Sys.time()
cat(sprintf("Feature counts: %0.2f minutes\n", difftime(featurecounts.end, featurecounts.start, units="mins")), file="benchmark.log", append=T)

featurecounts_rrna.start <- Sys.time()
if(file.exists(rrna_file)) {
    rrna_feature_counts = featureCounts(file.path(bam_folder, bam_file),
                                        annot.ext              = rrna_file,
                                        isGTFAnnotationFile    = TRUE,
                                        countMultiMappingReads = TRUE,
                                        fraction               = TRUE,
                                        strandSpecific         = stranded,
                                        isPairedEnd            = TRUE,
                                        nthreads               = n_threads)
} 
featurecounts_rrna.end <- Sys.time()
cat(sprintf("rRNA counts: %0.2f minutes\n", difftime(featurecounts_rrna.end, featurecounts_rrna.start, units="mins")), file="benchmark.log", append=T)


#########################################################################################################
### mapping counts

plotting.start <- Sys.time()
gene_count_stats <- gene_feature_counts$stat
colnames(gene_feature_counts$counts) = gsub(".bam", "", colnames(gene_feature_counts$counts))

rownames(gene_count_stats) <- gene_count_stats$Status
gene_count_stats <- gene_count_stats[,seq(2, ncol(gene_count_stats))]
colnames(gene_count_stats) <- gsub(".bam", "", colnames(gene_count_stats))

transformed_stats <- melt(t(gene_count_stats))

colnames(transformed_stats) <- c("Sample", "Group", "Alignments")

transformed_stats$Group <- factor(transformed_stats$Group,
                                 levels = rev(levels(transformed_stats$Group)[order(levels(transformed_stats$Group))]))
transformed_stats$Sample <- factor(transformed_stats$Sample,
                                  levels = rev(levels(transformed_stats$Sample)[order(as.character(transformed_stats$Sample))]))

transformed_stats <- transformed_stats[transformed_stats$Alignments > 0,]
transformed_stats$Reference = "Genes"


#########################################################################################################
# rRNA counts

if(file.exists(rrna_file)) {
    rrna_count_stats <- rrna_feature_counts$stat
    colnames(rrna_count_stats) <- gsub(".bam", "", colnames(rrna_count_stats))
    
    rrna_counts = as.data.frame(t(rrna_count_stats[rrna_count_stats$Status == "Assigned",2:ncol(rrna_count_stats)]))
    colnames(rrna_counts) = "Alignments"
    
    rrna_counts$Sample = rownames(rrna_counts)
    rrna_counts$Group = "rRNA"
    rrna_counts$Reference = "rRNA"
    
    d = rbind(transformed_stats, rrna_counts)
} else {
    d = transformed_stats
}

rrna_contamination_plot = ggplot(data = d, aes(x = Sample, y = Alignments)) +
  geom_col(aes(fill = Group), width = 0.7) +
  theme_bw() + facet_wrap(~Reference) + 
  ylab("Number of Alignments") + xlab("Samples") +
  ggtitle("Read Mapping Numbers") +
  theme(text = element_text(size = 18)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0))
print(rrna_contamination_plot)

write.table(d, "results/counts/mapping_stats.tsv", quote=F, sep="\t", row.names=F)

#########################################################################################################
# reads per gene

# get the reads per gene counts
gene_count_matrix = gene_feature_counts$counts
# set the sample names as column names
colnames(gene_count_matrix) = gsub(".bam", "", colnames(gene_count_matrix))

# group and count genes in classes of specific read counts
read_count_classes = data.frame("no_reads"=colSums(gene_count_matrix == 0),
                                "at_least_1_read"=colSums(gene_count_matrix >= 1 & gene_count_matrix < 10),
                                "at_least_10_reads"=colSums(gene_count_matrix >= 10 & gene_count_matrix < 100),
                                "at_least_100_reads"=colSums(gene_count_matrix >= 100 & gene_count_matrix < 1000),
                                "at_least_1000_reads"=colSums(gene_count_matrix >= 1000))

# setup data.frame for plotting
read_count_classes$Sample = rownames(read_count_classes)
melt_rcc = melt(read_count_classes)
melt_rcc$variable = as.character(melt_rcc$variable)

# sort the gene groups
melt_rcc$variable = factor(melt_rcc$variable, levels=c("no_reads",
                                                       "at_least_1_read",
                                                       "at_least_10_reads",
                                                       "at_least_100_reads",
                                                       "at_least_1000_reads"))

# plot the data as bar plot
gene_coverage_plot <- ggplot(melt_rcc, aes(x=Sample, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  ylab("Number of Genes") + xlab("Samples") +
  ggtitle("Number of Reads per Gene") +
  guides(fill=guide_legend(title="Number of assigned reads")) +
  theme(text = element_text(size = 18)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0))
print(gene_coverage_plot)

write.table(melt_rcc, "results/counts/gene_capture_stats.tsv", quote=F, sep="\t", row.names=F)

#########################################################################################################
# normalization

Length_kb = gene_feature_counts$annotation$Length / 1000
RPK = gene_feature_counts$counts / Length_kb
scaling_factors = colSums(RPK) / 1e6
TPM = RPK / scaling_factors

write.table(gene_feature_counts$counts, file=file.path(counts_folder, "feature.counts.tsv"), quote=FALSE, sep='\t')
write.table(TPM, file=file.path(counts_folder, "TPM.normalized.tsv"), quote=FALSE, sep='\t')

print(pheatmap(cor(log2(TPM+1)), fontsize = 18, main="Sample Correlation Heatmap"))
plotting.end <- Sys.time()

write.table(cor(log2(TPM+1)), "results/counts/sample_correlation.tsv", quote=F, sep="\t", row.names=F)

dev.off()

cat(sprintf("Plotting: %0.2f minutes\n", difftime(plotting.end, plotting.start, units="mins")), file="benchmark.log", append=T)




