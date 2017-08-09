############################################################
#### Script to perform differential expression analysis ####
############################################################

setwd("/Users/eling01/Google Drive/")

library(scran)

# Load normalized data
exp.data <- read.table('Platy/6th_approach/Data/Norm/norm_data.txt', sep = '\t', header = TRUE)

# Load cluster identities
clusters <- read.table("Platy/6th_approach/Results/stable_clusters_merged.txt", sep = "\t")

# Find marker genes for undifferentiated versus differentiated cells
clusters.1 <- data.frame(row.names = rownames(clusters), col = ifelse(grepl("Undiff", clusters$clusters), 1, 2), clusters = ifelse(grepl("Undiff", clusters$clusters), "Undiff", "Diff"))
markers <- findMarkers(as.matrix(log2(exp.data[,rownames(clusters)] + 1)), clusters.1$clusters)

# Select genes with FDR < 0.1 and logFC > 0
markers.spec <- lapply(markers, function(n){
  cur_n <- n[as.numeric(n[,3]) < 0.1 & as.numeric(n[,4]) > 0,]
  cur_n <- cur_n[order(as.numeric(cur_n$FDR), decreasing = FALSE),]
  cur_n
})

write.csv(as.data.frame(markers.spec$Diff), "Platy/6th_approach/Results/Marker_genes/Diff_markers.csv")
write.csv(as.data.frame(markers.spec$Undiff), "Platy/6th_approach/Results/Marker_genes/Undiff_markers.csv")

# Find marker genes for each differentiated cell group
clusters.2 <- clusters[!grepl("Undiff", clusters$clusters),]
markers <- findMarkers(as.matrix(log2(exp.data[,rownames(clusters.2)] + 1)), as.character(clusters.2$clusters))

# Select genes with FDR < 0.1 and logFC > 0
markers.spec <- lapply(markers, function(n){
  cur_n <- n[as.numeric(n[,3]) < 0.1 & as.numeric(n[,4]) > 0 & as.numeric(n[,5]) > 0 & as.numeric(n[,6]) > 0 & as.numeric(n[,7]) > 0,]
  cur_n <- cur_n[order(as.numeric(cur_n$FDR), decreasing = FALSE),]
  cur_n
})

write.csv(as.data.frame(markers.spec$ANS), "Platy/6th_approach/Results/Marker_genes/ANS_markers.csv")
write.csv(as.data.frame(markers.spec$Cilia), "Platy/6th_approach/Results/Marker_genes/Cilia_markers.csv")
write.csv(as.data.frame(markers.spec$DiffTrunk), "Platy/6th_approach/Results/Marker_genes/DiffTrunk_markers.csv")
write.csv(as.data.frame(markers.spec$Gut), "Platy/6th_approach/Results/Marker_genes/Gut_markers.csv")
write.csv(as.data.frame(markers.spec$Muscle), "Platy/6th_approach/Results/Marker_genes/Muscle_markers.csv")
