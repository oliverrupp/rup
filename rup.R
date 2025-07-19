library(ggplot2)
library(fastqcr)
library(Rsubread)
library(reshape2)
library(tidyr)
library(pheatmap)
library(Rfastp)
library(getopt)

theme_set(theme_bw())

spec = matrix(c(
  'datafolder' , 'd', 1, "character",
  'threads', 't', 1, "integer",
  'help', 'h', 0, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

n_threads     <- 1
source_folder <- "."

if (!is.null(opt$threads)) n_threads <- opt$threads
if (!is.null(opt$datafolder)) source_folder <- opt$datafolder
message("")

setwd(source_folder)

reference_folder      <- file.path(source_folder, "reference")
read_file_folder      <- file.path(source_folder, "reads")
results_folder        <- file.path(source_folder, "results")

genome_fasta_file     <- file.path(reference_folder, "genome.fa")
annotation_file       <- file.path(reference_folder, "annotation.gtf")
rrna_file             <- file.path(reference_folder, "rRNA.gtf")
samples_file          <- file.path(reference_folder, "samples.txt")

fastqc_folder         <- file.path(results_folder,   "fastqc")
trimmed_fastqc_folder <- file.path(results_folder,   "trimmed_fastqc")
trimmed_read_folder   <- file.path(results_folder,   "trimmed")
bam_folder            <- file.path(results_folder,   "bam")

sample_prefixes       <- gsub("_1.fq.gz", "", (list.files(read_file_folder, pattern = "*_1.fq.gz")))

if(!file.exists(annotation_file)) {
    message("annotation file missing [annotation.gtf]")
    quit(status=1)
}


for(folder in c(results_folder, fastqc_folder, trimmed_fastqc_folder, trimmed_read_folder, bam_folder)) {
  if(!dir.exists(folder)) {
    dir.create(folder)
  }
}

pdf("RNAseq_QC.pdf", w=18, h=12)


#########################################################################################################
#
# read based statistics
#
#########################################################################################################

#########################################################################################################
#### fastqc raw data

message("running fastqc ...")
### fastqc(fq.dir=read_file_folder, qc.dir=fastqc_folder, threads=n_threads)

qc <- qc_aggregate(fastqc_folder, progress=F)
qc$tot.seq <- as.numeric(qc$tot.seq)
qc$sample = gsub(".f(ast)?q.gz", "" , qc$sample)


#########################################################################################################
#### trimming

for(prefix in sample_prefixes) {
  outputPrefix <- file.path(trimmed_read_folder, prefix)

  if(!file.exists(paste(outputPrefix,"_R1.fastq.gz", sep=""))) {
     message(paste(outputPrefix,"_R1.fastq.gz", sep=""))
    fastp_stats <- rfastp(read1 = file.path(read_file_folder, paste(prefix, "_1.fq.gz", sep="")),
                          read2 = file.path(read_file_folder, paste(prefix, "_2.fq.gz", sep="")),
			  minReadLength = 25,
                          outputFastq = outputPrefix, thread = n_threads)
  }
}



#########################################################################################################
#### fastqc trimmed data
### fastqc(fq.dir=trimmed_read_folder, qc.dir=trimmed_fastqc_folder, threads=n_threads)

qc_trimmed <- qc_aggregate(trimmed_fastqc_folder, progress=F)
qc_trimmed$tot.seq <- as.numeric(qc_trimmed$tot.seq)
qc_trimmed$sample = gsub("_R([12])", "_\\1" , qc_trimmed$sample)

qc$Trimming = "raw"
qc_trimmed$Trimming = "trimmed"

qc_all = rbind(qc, qc_trimmed)

message(head(qc_all))

read_number_plot <- ggplot(qc_all, aes(x=sample, y=tot.seq, fill=Trimming)) + 
  geom_bar(stat="identity", position = "dodge") + 
  xlab("Samples") + ylab("Number of Reads") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  theme(text = element_text(size = 18)) + 
  ggtitle("Number of reads before and afer trimming")
print(read_number_plot)
 

#########################################################################################################
#
# mapping based statistics
#
#########################################################################################################

#########################################################################################################
### read alignment

if(!file.exists(file.path(reference_folder,"subread.index.reads"))) {
  buildindex(file.path(reference_folder,"subread.index"), genome_fasta_file)
}

for(prefix in sample_prefixes) {
  message(prefix)

  output_bam = file.path(bam_folder, paste(prefix, ".bam", sep=""))

  if(!file.exists(output_bam)) {
    mapping_stats <- align(index=file.path(reference_folder,"subread.index"),
                           readfile1=file.path(trimmed_read_folder, paste(prefix, "_R1.fastq.gz", sep="")),
                           readfile2=file.path(trimmed_read_folder, paste(prefix, "_R2.fastq.gz", sep="")),
                           output_file=output_bam, type=0,
                           minFragLength=0, maxFragLength=300, PE_orientation="fr",
                           nthreads=n_threads,
                           useAnnotation=TRUE,
                           annot.ext=annotation_file,
                           isGTF=TRUE,
                           nBestLocations=2)
  }
}





#########################################################################################################
### read counting

bam_file <- list.files(bam_folder, pattern = "*.bam$")

gene_feature_counts = featureCounts(file.path(bam_folder, bam_file),
                                    annot.ext              = annotation_file,
                                    isGTFAnnotationFile    = TRUE,
                                    countMultiMappingReads = FALSE,
                                    strandSpecific         = 0,
                                    isPairedEnd            = TRUE,
                                    nthreads               = n_threads)


if(file.exists(rrna_file)) {
    rrna_feature_counts = featureCounts(file.path(bam_folder, bam_file),
                                        annot.ext              = rrna_file,
                                        isGTFAnnotationFile    = TRUE,
                                        countMultiMappingReads = TRUE,
                                        fraction               = TRUE,
                                        strandSpecific         = 0,
                                        isPairedEnd            = TRUE,
                                        nthreads               = n_threads)
} 



#########################################################################################################
### mapping counts

gene_count_stats <- gene_feature_counts$stat

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
 

#########################################################################################################
# reads per gene

gene_count_matrix = gene_feature_counts$counts
colnames(gene_count_matrix) = gsub(".bam", "", colnames(gene_count_matrix))

df = data.frame("G0"=colSums(gene_count_matrix == 0),
                "G1"=colSums(gene_count_matrix >= 1 & gene_count_matrix < 10),
                "G10"=colSums(gene_count_matrix >= 10 & gene_count_matrix < 100),
                "G100"=colSums(gene_count_matrix >= 100 & gene_count_matrix < 1000),
                "G1000"=colSums(gene_count_matrix >= 1000))

df$Sample = rownames(df)

dfm = melt(df)

dfm$variable = as.character(dfm$variable)

dfm[dfm$variable == "G0",]$variable = "No reads"
dfm[dfm$variable == "G1",]$variable = "1 to 10 reads"
dfm[dfm$variable == "G10",]$variable = "10 to 100 reads"
dfm[dfm$variable == "G100",]$variable = "100 to 1000 reads"
dfm[dfm$variable == "G1000",]$variable = "more than 1000 reads"

dfm$variable = factor(dfm$variable, levels=c("No reads", "1 to 10 reads", "10 to 100 reads", "100 to 1000 reads", "more than 1000 reads"))

gene_coverage_plot <- ggplot(dfm, aes(x=Sample, y=value, fill=variable)) + 
  geom_bar(stat="identity") +
  ylab("Number of Genes") + xlab("Samples") +
  ggtitle("Number of Reads per Gene") +
  guides(fill=guide_legend(title="Number of assigned reads")) + 
  theme(text = element_text(size = 18)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0))
print(gene_coverage_plot)

#########################################################################################################
# normalization

Length_kb = gene_feature_counts$annotation$Length / 1000
RPK = gene_feature_counts$counts / Length_kb
colnames(RPK) = gsub(".bam", "", colnames(RPK))
scaling_factors = colSums(RPK) / 1e6
TPM = RPK / scaling_factors

print(pheatmap(cor(log2(TPM+1)), fontsize = 18, main="Sample Correlation Heatmap"))

dev.off()


