library(ggplot2)
library(fastqcr)
library(Rsubread)
library(reshape2)
library(tidyr)
library(pheatmap)
library(Rfastp)

theme_set(theme_bw())

n_threads             <- 32 

source_folder         <- "/var/scratch/orupp/RNAseqQC"

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
fastqc(fq.dir=read_file_folder, qc.dir=fastqc_folder, threads=n_threads)

qc <- qc_aggregate(fastqc_folder)
qc$tot.seq <- as.numeric(qc$tot.seq)


# read_number_plot <- ggplot(qc[qc$module == "Basic Statistics",], aes(x=sample, y=tot.seq)) + geom_bar(stat="identity") + ggtitle("raw read number")
# read_pct_dup_plot <- ggplot(qc[qc$module == "Basic Statistics",], aes(x=sample, y=pct.dup)) + geom_bar(stat="identity") + ggtitle("raw read duplication")

# print(read_number_plot)
# print(read_pct_dup_plot)



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
fastqc(fq.dir=trimmed_read_folder, qc.dir=trimmed_fastqc_folder, threads=n_threads)

qc_trimmed <- qc_aggregate(trimmed_fastqc_folder)
qc_trimmed$tot.seq <- as.numeric(qc_trimmed$tot.seq)

# read_number_plot <- ggplot(qc_trimmed[qc_trimmed$module == "Basic Statistics",], aes(x=sample, y=tot.seq)) + geom_bar(stat="identity") + ggtitle("trimmed reaad number")
# read_pct_dup_plot <- ggplot(qc_trimmed[qc_trimmed$module == "Basic Statistics",], aes(x=sample, y=pct.dup)) + geom_bar(stat="identity") + ggtitle("trimmed read duplication")

qc$reads = "raw"
qc_trimmed$reads = "trimmed"

qc_all = rbind(qc, qc_trimmed)

read_number_plot <- ggplot(qc_all[qc_all$module == "Basic Statistics",], aes(x=sample, y=tot.seq, fill=reads)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle("read number")

print(read_number_plot)


#########################################################################################################
#
# mapping based statistics
#
#########################################################################################################

#########################################################################################################
### read alignment

buildindex(file.path(reference_folder,"subread.index"), genome_fasta_file)

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

gene_counts = featureCounts(file.path(bam_folder, bam_file),
                            annot.ext              = annotation_file,
                            isGTFAnnotationFile    = TRUE,
                            countMultiMappingReads = FALSE,
                            strandSpecific         = 0,
                            isPairedEnd            = TRUE,
                            nthreads               = n_threads)


rrna_counts = featureCounts(file.path(bam_folder, bam_file),
                            annot.ext              = rrna_file,
                            isGTFAnnotationFile    = TRUE,
                            countMultiMappingReads = TRUE,
                            fraction               = TRUE,
                            strandSpecific         = 0,
                            isPairedEnd            = TRUE,
                            nthreads               = n_threads)

#rrna_counts_u = featureCounts(file.path(bam_folder, bam_file),
#                            annot.ext              = rrna_file,
#                            isGTFAnnotationFile    = TRUE,
#                            countMultiMappingReads = FALSE,
#                            fraction               = FALSE,
#                            strandSpecific         = 0,
#                            isPairedEnd            = TRUE,
#                            nthreads               = n_threads)




#########################################################################################################
### mapping counts

gene_count_stats <- gene_counts$stat
rrna_count_stats <- rrna_counts$stat

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

#mapping_stats_plot <- ggplot(data = transformed_stats, aes(x = Sample, y = Alignments)) +
#  geom_col(aes(fill = Group), width = 0.7)+
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = -22.5, hjust = 0, size=12))
#print(mapping_stats_plot)



#########################################################################################################
# rRNA counts

rrna_counts = as.data.frame(t(rrna_count_stats[rrna_count_stats$Status == "Assigned",2:ncol(rrna_count_stats)]))
colnames(rrna_counts) = "Alignments"

rrna_counts$Sample = rownames(rrna_counts)
rrna_counts$Group = "rRNA"

# ggplot(rrna_counts, aes(y = Alignments, x = Sample)) + geom_bar(stat = "identity")


d = rbind(transformed_stats, rrna_counts)

rrna_contamination_plot = ggplot(data = d, aes(x = Sample, y = Alignments)) +
  geom_col(aes(fill = Group), width = 0.7)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = -22.5, hjust = 0, size=12))
print(rrna_contamination_plot)



#########################################################################################################
# reads per gene

gene_count_matrix = gene_counts$counts

df = data.frame("0"=colSums(gene_count_matrix == 0),
                "1"=colSums(gene_count_matrix >= 1 & gene_count_matrix < 10),
                "10"=colSums(gene_count_matrix >= 10 & gene_count_matrix < 100),
                "100"=colSums(gene_count_matrix >= 100 & gene_count_matrix < 1000),
                "1000"=colSums(gene_count_matrix >= 1000))

df$Sample = rownames(df)

gene_coverage_plot <- ggplot(melt(df), aes(x=Sample, y=value, fill=variable)) + geom_bar(stat="identity")
print(gene_coverage_plot)

#########################################################################################################
# normalization

Length_kb = gene_counts$annotation$Length / 1000
RPK = gene_counts$counts / Length_kb
scaling_factors = colSums(RPK) / 1e6
TPM = RPK / scaling_factors


#top_TPM = TPM[rowSums(TPM) > 10 & sum(TPM > 5) > 5,]

print(pheatmap(cor(log2(TPM+1))))


# tpm_centered <- t(TPM-rowMeans(TPM))
# tpm_svd <- svd(tpm_centered)
# plot(tpm_svd$u[,1], tpm_svd$u[,2])
# tpm_prcomp <- prcomp(tpm_centered)
# plot(tpm_prcomp$x[,1], tpm_prcomp$x[,2])

# head(tpm_prcomp$x[,1])

# plot_df <- data.frame(PC1 = tpm_prcomp$x[,1], PC2 = tpm_prcomp$x[,2], Samples = rownames(tpm_prcomp$x))
# pca_plot <- ggplot(plot_df, aes(x = PC1, y = PC2, col = Samples)) + geom_point()


dev.off()


