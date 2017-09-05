###################################################
#### Script to recreate extended data figure 2 ####
###################################################

library(scran)
library(Rtsne)
library(plot3D)

# Load normalized data
exp.data <- read.table("/Users/eling01/Google Drive/Platy/6th_approach/Data/Norm/norm_data.txt", sep = "\t")

# Load cluster identities
clusters <- read.table("/Users/eling01/Google Drive/Platy/6th_approach/Results/stable_clusters.txt", sep = "\t", stringsAsFactors = FALSE)

# Find group specifc markers for all 7 groups
markers <- findMarkers(log(as.matrix(exp.data[,rownames(clusters)] + 1)), as.character(clusters$clusters))

# Load highly variable genes
HVG <- read.table("/Users/eling01/Google Drive/Platy/6th_approach/Results/HVG.txt", sep = "\t")

# EDF2 a Visualize groups in tSNE
set.seed(3)
tsne <- Rtsne(log10(t(exp.data[as.character(HVG$GeneNames[HVG$HVG == TRUE]),rownames(clusters)] + 1)), perplexity = 30)
plot(tsne$Y[,1], tsne$Y[,2], pch = 16, col = c("red", "yellow", "dark grey", "blue", "green", "brown", "purple")[as.factor(clusters$clusters)])

# EDF2 b-h Visualize closeness of groups via logFC
pdf("/Users/eling01/Google Drive/Platy/6th_approach/Figures/Supplements/EDF2bh.pdf")
plot(hclust(dist(t(markers$ANS[markers$ANS$FDR < 0.1,4:9]))), hang = -1)
plot(hclust(dist(t(markers$Cilia[markers$Cilia$FDR < 0.1,4:9]))), hang = -1)
plot(hclust(dist(t(markers$DiffTrunk[markers$DiffTrunk$FDR < 0.1,4:9]))), hang = -1)
plot(hclust(dist(t(markers$Gut[markers$Gut$FDR < 0.1,4:9]))), hang = -1)
plot(hclust(dist(t(markers$Muscle[markers$Muscle$FDR < 0.1,4:9]))), hang = -1)
plot(hclust(dist(t(markers$UndiffHead[markers$UndiffHead$FDR < 0.1,4:9]))), hang = -1)
plot(hclust(dist(t(markers$UndiffTrunk[markers$UndiffTrunk$FDR < 0.1,4:9]))), hang = -1)
dev.off()

# Two groups - UndiffHead and UndiffTrunk are consistently closest to each other
# Merge them to one group

# EDF2 i Visualize PCNA expression as marker gene for differentiating cells
gene <- "PCNA"
library(vioplot)
vioplot(log10(exp.data[gene,rownames(clusters)[clusters$clusters == "ANS"]] + 1),
        log10(exp.data[gene,rownames(clusters)[clusters$clusters == "Cilia"]] + 1),
        log10(exp.data[gene,rownames(clusters)[clusters$clusters == "Muscle"]] + 1),
        log10(exp.data[gene,rownames(clusters)[clusters$clusters == "Gut"]] + 1),
        log10(exp.data[gene,rownames(clusters)[clusters$clusters == "DiffTrunk"]] + 1),
        log10(exp.data[gene,rownames(clusters)[clusters$clusters == "UndiffHead"]] + 1),
        log10(exp.data[gene,rownames(clusters)[clusters$clusters == "UndiffTrunk"]] + 1))
