######################################
#### Script to reproduce Figure 1 #### 
######################################

library(Rtsne)
library(pheatmap)

# Load normalized data
exp.data <- read.table("/Users/nils/Google Drive/Platy/6th_approach/Data/Norm/norm_data.txt", sep = "\t")

# Load cluster identities
clusters <- read.table("/Users/nils/Google Drive/Platy/6th_approach/Results/stable_clusters_merged.txt", sep = "\t")

# Load highly variable genes
HVG <- read.table("/Users/nils/Google Drive/Platy/6th_approach/Results/HVG.txt", sep = "\t")

# Calculate tSNE on highly variable genes
set.seed(3)
tsne <- Rtsne(log10(t(exp.data[as.character(HVG$GeneNames[HVG$HVG == TRUE]),rownames(clusters)] + 1)), perplexity = 30)

# Plot tSNE
plot(tsne$Y[,1], tsne$Y[,2], pch = 21, bg = c("red", "yellow", "dark grey", "blue", "green", "white")[clusters$clusters])

# Visualise the molecular features (specifc genes) of each cluster as determined by sparse clustering
load("/Users/nils/Google Drive/Platy/6th_approach/Results/Iter_clust.RData")

genes <- names(groups$group7[[1]]$ws)[groups$group7[[1]]$ws > 0.03]
pheatmap(log10(exp.data[genes,c(rownames(clusters)[clusters$clusters == "ANS"], rownames(clusters)[clusters$clusters == "Gut"], rownames(clusters)[clusters$clusters == "Muscle"], rownames(clusters)[clusters$clusters == "Cilia"], rownames(clusters)[clusters$clusters == "DiffTrunk"])] + 1), 
         cluster_cols = FALSE, cellwidth = 2, cellheight = 2, 
         col = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100), 
         show_rownames = FALSE, show_colnames = FALSE)     
