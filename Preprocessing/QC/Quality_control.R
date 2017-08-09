###################################################################
#### Script to perform quality control on 48hpf data ##############
###################################################################

library(readxl)
library(RColorBrewer)
library(ggplot2)

setwd("/Users/eling01/Google Drive/")

# Read in metadata on cells
QC <- as.data.frame(read_xlsx("Platy/6th_approach/Data/QC/metadata.xlsx"))

# Read in the transcript counts
raw.data <- read.table("Platy/6th_approach/Data/Raw/170502_annotGenes.txt", sep = "\t", header = TRUE)

# Select annotated cells
raw.data <- raw.data[,colnames(raw.data) %in% as.character(QC$SampleID)]
QC <- QC[as.character(QC$SampleID) %in% colnames(raw.data),]

# Fig S1A - Plot number of mapped reads, color by empty well, single of multiple
plot(log10(colSums(raw.data)), pch = 16, col = ifelse(as.numeric(QC$nr.of.cells) == 1 & is.na(QC$note), "goldenrod", ifelse(as.numeric(QC$nr.of.cells) > 1 & is.na(QC$note), "darkgreen", ifelse(as.numeric(QC$nr.of.cells) == 0 & is.na(QC$note), "darkblue", "black"))))

# Select single cells 
single <- as.character(QC$SampleID[QC$nr.of.cells == 1 & is.na(QC$note)])
raw.data <- raw.data[,single]

# Fig S1B - ERCC reads vs transcriptomic reads
plot(log10(colSums(raw.data[!grepl("ERCC-0", rownames(raw.data)),])), 
     log10(colSums(raw.data[grepl("ERCC-0", rownames(raw.data)),])), pch = 16,
     col = ifelse(log10(colSums(raw.data[!grepl("ERCC-0", rownames(raw.data)),])) < 5 | log10(colSums(raw.data[grepl("ERCC-0", rownames(raw.data)),])) < 3, "red", "black"))

raw.data <- raw.data[log10(colSums(raw.data[grepl("ERCC-0", rownames(raw.data)),])) > 3 & log10(colSums(raw.data[!grepl("ERCC-0", rownames(raw.data)),])) > 5]

# Fig S1C - Percentage of reads mapping to transcriptome
plot(colSums(raw.data[!grepl("ERCC-0", rownames(raw.data)),])/colSums(raw.data), col = ifelse(colSums(raw.data[!grepl("ERCC-", rownames(raw.data)),])/colSums(raw.data) < 0.6 , "red", "black"), pch = 16)
raw.data <- raw.data[which(colSums(raw.data[!grepl("ERCC-0", rownames(raw.data)),])/colSums(raw.data) > 0.60)]

# Fig S1D - Number of genes expressed in each cell
plot(apply(raw.data, 2, function(n){length(which(n > 0))}), pch = 16, col = ifelse(apply(raw.data, 2, function(n){length(which(n > 0))}) < 2000, "red", "black"), ylab = "No. genes detected")
raw.data <- raw.data[,apply(raw.data, 2, function(n){length(which(n > 0))}) > 2000]

# Fig S1E - Remove lowly and highly expressed genes and ERCCs that are not expressed
ERCC <- raw.data[grepl("ERCC-0", rownames(raw.data)),]
ERCC <- ERCC[rowSums(ERCC)>0,]
counts <- raw.data[!grepl("ERCC-0", rownames(raw.data)),]
plot(log10(rowSums(counts) + 1), pch = 16, cex = 0.5, col = ifelse(log10(rowSums(counts) + 1) < 1 | log10(rowSums(counts) + 1) > 6, "purple", "black"))
counts <- counts[log10(rowSums(counts) + 1) > 1 & log10(rowSums(counts) + 1) < 6,]

raw.data <- rbind(counts, ERCC)

write.table(as.data.frame(raw.data), "/Users/eling01/Google Drive/Platy/6th_approach/Data/Raw/QC_rawData.tsv", sep = "\t")

# Fig S1G - Stats on read mapping
# Read in the initial raw transcript counts
raw.data.init <- read.table("Platy/6th_approach/Data/Raw/CN61_genes.raw.htseq1.tsv", sep = "\t", header = TRUE)
raw.data.init <- raw.data.init[,colnames(raw.data)]

htseq.stats <- read.table("Platy/6th_approach/Data/all_HTSeqStats.txt", sep = "\t", header = TRUE)
htseq.stats <- htseq.stats[,colnames(raw.data)]

# Collect relevant stats
df <- data.frame(row.names = colnames(htseq.stats), 
                 total_reads = as.numeric(colSums(rbind(raw.data.init, htseq.stats))), 
                 mapped_reads = as.numeric(colSums(raw.data.init)), 
                 ambig = as.numeric(htseq.stats[2,]), 
                 multi = as.numeric(htseq.stats[5,]))

# Ratio to all reads
df[,2:4] <- df[,2:4]/df[,1]

barplot(as.matrix(t(df[,2:4])), col = c("black", "white", "orange"))
