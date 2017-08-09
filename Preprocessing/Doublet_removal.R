###################################
#### Script to remove doublets ####
###################################

# Load normalized data
norm.data <- read.table("/Users/nils/Google Drive/Platy/6th_approach/Data/Norm/norm_data.txt", header = TRUE, sep = "\t")

# Load cluster identities
clusters <- read.table("/Users/nils/Google Drive/Platy/6th_approach/Results/stable_clusters.txt", header = TRUE, sep = "\t")

# Read in marker genes for each cluster
ANS.genes <- read.csv("/Users/nils/Google Drive/Platy/6th_approach/Results/Marker_genes/ANS_markers.csv")
Cilia.genes <- read.csv("/Users/nils/Google Drive/Platy/6th_approach/Results/Marker_genes/Cilia_markers.csv")
Gut.genes <- read.csv("/Users/nils/Google Drive/Platy/6th_approach/Results/Marker_genes/Gut_markers.csv")
Muscle.genes <- read.csv("/Users/nils/Google Drive/Platy/6th_approach/Results/Marker_genes/Muscle_markers.csv")
Trunk.genes <- read.csv("/Users/nils/Google Drive/Platy/6th_approach/Results/Marker_genes/DiffTrunk_markers.csv")

# Doublets should be characterized by a high expression of marker genes of a different groups
# This can be visualized by boxplots
# Remove doublets from differentiated cell groups
clusters <- clusters[!grepl("Undiff", clusters$clusters),]
norm.data <- norm.data[,rownames(clusters)]
barplot(apply(norm.data[as.character(ANS.genes$Gene),rownames(clusters)[!(clusters$clusters == "ANS")]], 2, function(n){length(which(n>0))})/apply(norm.data[,rownames(clusters)[!(clusters$clusters == "ANS")]], 2, function(n){length(which(n>0))}))
ANS.doublets <- rownames(clusters)[!(clusters$clusters == "ANS")][apply(norm.data[as.character(ANS.genes$Gene),rownames(clusters)[!(clusters$clusters == "ANS")]], 2, function(n){length(which(n>0))})/apply(norm.data[,rownames(clusters)[!(clusters$clusters == "ANS")]], 2, function(n){length(which(n>0))}) > 0.04]

barplot(apply(norm.data[as.character(Cilia.genes$Gene),rownames(clusters)[!(clusters$clusters == "Cilia")]], 2, function(n){length(which(n>0))})/apply(norm.data[,rownames(clusters)[!(clusters$clusters == "Cilia")]], 2, function(n){length(which(n>0))}))

barplot(apply(norm.data[as.character(Muscle.genes$Gene),rownames(clusters)[!(clusters$clusters == "Muscle")]], 2, function(n){length(which(n>0))})/apply(norm.data[,rownames(clusters)[!(clusters$clusters == "Muscle")]], 2, function(n){length(which(n>0))}))

barplot(apply(norm.data[as.character(Gut.genes$Gene),rownames(clusters)[!(clusters$clusters == "Gut")]], 2, function(n){length(which(n>0))})/apply(norm.data[,rownames(clusters)[!(clusters$clusters == "Gut")]], 2, function(n){length(which(n>0))}))

barplot(apply(norm.data[as.character(Trunk.genes$Gene),rownames(clusters)[!(clusters$clusters == "DiffTrunk")]], 2, function(n){length(which(n>0))})/apply(norm.data[,rownames(clusters)[!(clusters$clusters == "DiffTrunk")]], 2, function(n){length(which(n>0))}))

# Remove the doublet
norm.data <- read.table("/Users/eling01/Google Drive/Platy/6th_approach/Data/Norm/norm_data.txt", header = TRUE, sep = "\t")
clusters <- read.table("/Users/eling01/Google Drive/Platy/6th_approach/Results/stable_clusters_170510.txt", header = TRUE, sep = "\t")

norm.data <- norm.data[,!grepl(ANS.doublets, colnames(norm.data))]
write.table(as.data.frame(norm.data), "/Users/eling01/Google Drive/Platy/6th_approach/Data/Norm/norm_data.txt", sep = "\t")
clusters <- clusters[!grepl(ANS.doublets, rownames(clusters)),]
write.table(as.data.frame(clusters), "/Users/eling01/Google Drive/Platy/6th_approach/Results/stable_clusters.txt", sep = "\t")
