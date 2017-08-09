####################################
#### Script to create Figure 5 ####
####################################

library(pheatmap)

setwd("/Users/eling01/Google Drive/")

# Load normalized data
exp.data <- read.table("Platy/6th_approach/Data/Norm/norm_data.txt", sep = "\t")
HVG <- read.table("Platy/6th_approach/Results/HVG.txt", sep = "\t", stringsAsFactors = FALSE)

# Load cluster identities
clusters <- read.table("Platy/6th_approach/Results/stable_clusters.txt", sep = "\t")

# Select genes for visualization
genes <- c("Syt1", "Rab3_complete", "7B2_neuropeptide_precursor_secretogranin-V", "Pdu_m-unc13_Contig13855", "RIMSBP_comp424246", "Syntaxin1a", "Syntaxin18__Loc_9481", "ERC1_ERC2_Cast", "Complexin_v2", "RIMS", "Synapsin3", "SNAP25",  "TBR1_MOUSE T-box brain protein 1", "Contig1452_Neuroglobin", "NPY_Primr", "Loc_6730_VAT1L", "cGMP-PDE_Loc_72981", "Sytalpha", "Syt17", "Syt12", "Syt4",  "PHC2",
           "g587-t1_cirri_antenna_parapodia_palpae_7TM3_GRM7_Contig11280", "NRT_DROME Neurotactin", "Pdu_similar-to-neurotrypsin", "Pdu_trscr_assembly_1|102268", "cand.HSPG2_Ischia", "Tbx20_Loc_14134", "uncx4", "Phox2")

# Load gene annotation
gene.names <- read.table(file = "/Users/eling01/Dropbox/Platy_neurons/Results/allgenes_abbr.txt", sep = "\t", header = TRUE, fill = TRUE, quote = "", stringsAsFactors = FALSE)

# Find gene names for annotation
annot.names <- gene.names[match(genes, gene.names[,1]),2]

# Match with Jekely annotation
GJ.genes <- read.table("Platy/6th_approach/Data/Gene_annotation/names_GJ_split.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
GJ.genes$Evalue <- as.numeric(sapply(sapply(as.character(GJ.genes$cut_blasthit), function(n){unlist(strsplit(n, split = "E-value: "))[2]}), function(x){unlist(strsplit(x, split = " "))[1]}))
GJ.genes$BLAST_hit <- as.character(sapply(sapply(as.character(GJ.genes$cut_blasthit), function(n){unlist(strsplit(n, split = "\\|"))[3]}), function(x){unlist(strsplit(x, split = " OS"))[1]}))

annot.names[is.na(annot.names)] <- gene.names[match(GJ.genes$jekely_ID[match(genes[which(is.na(annot.names))], GJ.genes$BLAST_hit)], gene.names[,1]),2]

# Figure 5B - Plot heatmap
pheatmap(log10(exp.data[genes,rownames(clusters)[-which(grepl("Undiff", clusters$clusters))]] + 1), 
         cluster_rows = FALSE, 
         annotation_col = data.frame(row.names = rownames(clusters)[-which(grepl("Undiff", clusters$clusters))], group = clusters$clusters[-which(grepl("Undiff", clusters$clusters))]), 
         clustering_distance_cols = as.dist((1- cor(log10(exp.data[as.character(HVG$GeneNames[HVG$HVG == TRUE]), rownames(clusters)[-which(grepl("Undiff", clusters$clusters))]] + 1), method = "spearman"))/2), 
         col = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100), cellwidth = 5, cellheight = 5, gaps_row = c(12, 22), fontsize = 5,
         labels_row = annot.names)
