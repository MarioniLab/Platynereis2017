######################################
#### Script to reproduce Figure 6 ####
######################################

setwd("/Users/eling01/Google Drive/")

library(pvclust)

# Load normalized data
exp.data <- read.table("Platy/6th_approach/Data/Norm/norm_data.txt", sep = "\t")
HVG <- read.table("Platy/6th_approach/Results/HVG.txt", sep = "\t", stringsAsFactors = FALSE)

# Load cluster identities
clusters <- read.table("Platy/6th_approach/Results/stable_clusters.txt", sep = "\t")

input <- exp.data[,rownames(clusters)[!grepl("Undiff", clusters$clusters)]]
input <- input[HVG$GeneNames[HVG$HVG == TRUE],]
dendo <- hclust(as.dist(sqrt((1-cor(log10(input + 1), method = "spearman"))/2)), method = "complete")

library(RColorBrewer)
labelColors = c("red", "green", "blue", "black", "yellow")
# cut dendrogram in 5 clusters
clusMember = clusters[clusters$clusters == "ANS" | clusters$clusters == "DiffTrunk" | clusters$clusters == "Cilia" | clusters$clusters == "Gut" | clusters$clusters == "Muscle",1] - 2
names(clusMember) <- colnames(input)
# function to get color labels
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}
# using dendrapply
clusDendro = dendrapply(as.dendrogram(dendo), colLab)
# make plot
plot(clusDendro, main = "Full Dendrogram")

# Use the pvclust package to calculate p-values
pv.clust <- pvclust(log10(input + 1), method.dist = function(x){as.dist(sqrt((1-cor(x, method = "spearman"))/2))}, method.hclust = "complete", nboot = 1000)
plot(pv.clust, hang=-1, cex = 0.3)

# Cell type tree by averaging across groups
tree.df <- data.frame(row.names = HVG$GeneNames[HVG$HVG == TRUE], 
                      ANS = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "ANS"]]), 
                      Cilia = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "Cilia"]]),
                      Muscle = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "Muscle"]]),
                      Gut = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "Gut"]]),
                      Trunk = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "DiffTrunk"]]))

# Visualize the dendrogram using spearman dissimilarity
dendo <- hclust(as.dist(sqrt((1-cor(log10(tree.df + 1), method = "spearman"))/2)), method = "complete")
plot(dendo)

# pvclust package to perform bootstrapping
library(pvclust)
pv.clust <- pvclust(log10(input + 1), method.dist = function(x){as.dist(sqrt((1-cor(x, method = "spearman"))/2))}, method.hclust = "complete", nboot = 1000)
plot(pv.clust, hang = -1)

# pvclust on the averaged tree
# Cell type tree
tree.df <- data.frame(row.names = HVG$GeneNames[HVG$HVG == TRUE], 
                      ANS = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "ANS"]]), 
                      Cilia = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "Cilia"]]),
                      Muscle = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "Muscle"]]),
                      Gut = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "Gut"]]),
                      Trunk = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "DiffTrunk"]]))

pv.clust <- pvclust(log10(tree.df + 1), method.dist = function(x){as.dist(sqrt((1-cor(x, method = "spearman"))/2))}, method.hclust = "complete", nboot = 1000)
plot(pv.clust, hang = -1)
