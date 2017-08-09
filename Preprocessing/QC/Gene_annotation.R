###############################################################
#### Script to match BLAST hits to gene names #################
###############################################################

# Read in the data - gene annotation and raw count files
gene_annotation <- read.table("/Users/eling01/Google Drive/Platy/6th_approach/Data/Gene_annotation/T90LvsUniref_reciprocal_only.txt", header = TRUE, fill = TRUE, sep = "\t")
input <- read.table("/Users/eling01/Google Drive/Platy/6th_approach/Data/Raw/CN61_genes.raw.htseq1.tsv", header=TRUE, sep="\t", row.names = 1)

# Spike-in reads
ERCC <- input[which(grepl("ERCC-", rownames(input))),]

# Biological genes
input <- input[-which(grepl("ERCC-", rownames(input))),]

# Separate between Arendt genes, Jekely genes and already annotated genes
input.Loc <- input[which(grepl("^Loc_[0-9]+_Tr_", rownames(input))),]
input.trscr <- input[which(grepl("Pdu_trscr_assembly_1", rownames(input))),]
input.annot <- input[-c(which(grepl("^Loc_[0-9]+_Tr_", rownames(input))), which(grepl("Pdu_trscr_assembly_1", rownames(input)))),]

# Order BLAST hits by E value
gene_annotation <- gene_annotation[order(gene_annotation$T90.to.Uniref_eval, decreasing = FALSE),]

# Pick first hit regrading transcriptome range - hit with highest E value
gene_annotation_unique <- gene_annotation[match(unique(as.character(gene_annotation$T90)), as.character(gene_annotation$T90)),]

# Pick first hit regarding BLAST result
gene_annotation_unique <- gene_annotation_unique[match(unique(as.character(gene_annotation_unique$T90.to.Uniref_Hit_description)), as.character(gene_annotation_unique$T90.to.Uniref_Hit_description)),]

# Replace rownames in Arendt genes with BLAST hit
m <- match(rownames(input.Loc), as.character(gene_annotation_unique$T90))
rownames(input.Loc)[!is.na(m)] <- as.character(gene_annotation_unique$T90.to.Uniref_Hit_description[m[!is.na(m)]])

# Replace Jekely genes by BLAST hits
GJ.genes <- read.table("/Users/eling01/Google Drive/Platy/6th_approach/Data/Gene_annotation/names_GJ_split.txt", sep = "\t", header = TRUE)
GJ.genes$Evalue <- as.numeric(sapply(sapply(as.character(GJ.genes$cut_blasthit), function(n){unlist(strsplit(n, split = "E-value: "))[2]}), function(x){unlist(strsplit(x, split = " "))[1]}))
GJ.genes$BLAST_hit <- as.character(sapply(sapply(as.character(GJ.genes$cut_blasthit), function(n){unlist(strsplit(n, split = "\\|"))[3]}), function(x){unlist(strsplit(x, split = " OS"))[1]}))

# Order BLAST hits by Evalues
GJ.genes <- GJ.genes[order(GJ.genes$Evalue, decreasing = FALSE),]

# Pick first hit regarding BLAST result
GJ.genes.unique <- GJ.genes[match(unique(as.character(GJ.genes$BLAST_hit)), as.character(GJ.genes$BLAST_hit)),]

# Replace rownames in Jekely genes with BLAST hit
m <- match(rownames(input.trscr), as.character(GJ.genes.unique$jekely_ID))
rownames(input.trscr)[!is.na(m)] <- as.character(GJ.genes.unique$BLAST_hit[m[!is.na(m)]])

# Reassemble data matrix
raw_counts <- rbind(input.annot, input.Loc, input.trscr, ERCC)

write.table(as.data.frame(raw_counts), "/Users/eling01/Google Drive/Platy/6th_approach/Data/Raw/170502_annotGenes.txt", sep = "\t", quote = FALSE)
