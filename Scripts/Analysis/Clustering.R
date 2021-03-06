################################################################
#### Script for cluster detection of 48hpf normalized data ####
###############################################################

setwd("/Users/eling01/Google Drive/")

library(Rtsne)
library(sparcl)
library(RColorBrewer)
library(plot3D)
library(pheatmap)
library(mclust)
library(dynamicTreeCut)
library(cluster)
library(clusteval)

# Read in data
load("Platy/6th_approach/Data/Norm/norm_data.RData")

# Estimate the best number of clusters
# Iterative clutstering using the sparcl package
clust.dat <- t(log10(exp.data[as.character(HVG$Table$GeneNames[HVG$Table$HVG == TRUE]),] + 1))
clust.dat.scaled <- scale(clust.dat, TRUE, TRUE)

groups <- list()
tunes <- list()

for(k in 2:10){
  tune <- KMeansSparseCluster.permute(x = clust.dat.scaled, K = k, wbounds = seq(2, 100, 2), nperms = 15)
  tunes[[paste("group", k, sep="")]] <- tune
  group <- KMeansSparseCluster(x = clust.dat.scaled, K = k, wbounds = tune$bestw)
  groups[[paste("group", k, sep="")]] <- group
}

# Visualize clutsering on tSNE
scatter2D(tsne$Y[,1], tsne$Y[,2], colvar = groups$group10[[1]]$Cs, pch = 16, col = brewer.pal(10, "Set3"))

load("Platy/6th_approach/Results/Iter_clust.RData")

# Plot the average within cluster sum of squares
ss.all <- list()

for(x in 1:length(groups)){
  ss <- 0
  
  cells <- log10(exp.data[as.character(HVG$Table$GeneNames[HVG$Table$HVG == TRUE]),] + 1)
  ### SSE distance for all cells
  clust.mean <- rowMeans(cells)
  r.dist <- apply(cells, 2, function(n){(n-clust.mean)^2})
  cur_ss <- sum(r.dist)
  
  ss.all[["group1"]] <- cur_ss
  
  for(i in 1:length(unique(groups[[x]][[1]]$Cs))){
    cells <- log10(exp.data[as.character(HVG$Table$GeneNames[HVG$Table$HVG == TRUE]),which(groups[[x]][[1]]$Cs == i)] + 1)
    ### SSE distance for clutsers
    clust.mean <- rowMeans(cells)
    r.dist <- apply(cells, 2, function(n){(n-clust.mean)^2})
    cur_ss <- sum(r.dist)
    ss <- ss + cur_ss
  }
  ss.all[[x+1]] <- ss/length(unique(groups[[x]][[1]]$Cs))
}

# Supplementary Figure S1k
plot(unlist(ss.all), pch =16, type = "l")
points(unlist(ss.all), pch =16)

# Show how clusters change across different Ks
mat <- matrix(unlist(lapply(groups, function(n){n[[1]]$Cs})), nrow = length(groups), ncol = ncol(exp.data), byrow = TRUE)
# Order based on clustering at K=7
mat <- mat[,order(groups$group7[[1]]$Cs, decreasing = TRUE)]
# Plot heatmap
pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, col = rainbow(10))

# Compare sparse clustering to other approaches that don't depend on K
# Dynamic tree clustering

dist.all <- as.dist(sqrt((1 - cor(log10(exp.data[as.character(HVG$Table$GeneNames[HVG$Table$HVG == TRUE]),] + 1), method = "spearman"))/2))
dendro <- hclust(dist.all, method = "ward.D2")
ct <- cutreeDynamic(dendro = dendro, distM = as.matrix(dist.all), minClusterSize = 5, deepSplit = 0)
plot(ct)

plot(tsne$Y[,1], tsne$Y[,2], pch = 16, col=brewer.pal(12, "Set3")[ct])

# Compare clusterings of each group with clusters of DynamicTreeCut

sim.jac <- vector(length=9)
sim.rand <- vector(length=9)

for(i in 1:9){
  sim.jac[i] =  cluster_similarity(labels1 = groups[[i]][[1]]$Cs, labels2 = ct, similarity = "jaccard")
  sim.rand[i] =  cluster_similarity(labels1 = groups[[i]][[1]]$Cs, labels2 = ct, similarity = "rand")
}

# Compare clustering to finite normal mixture modelling
X <- t(log10(exp.data[as.character(HVG$Table$GeneNames[HVG$Table$HVG == TRUE]),] + 1))

BIC = mclustBIC(X)

summary(BIC)

mod1 = Mclust(X, x=BIC)

# Compare clustering between finite mixture and sparse clustering

sim.jac.mclust <- vector(length=9)
sim.rand.mclust <- vector(length=9)

for(i in 1:9){
  sim.jac.mclust[i] =  cluster_similarity(labels1 = groups[[i]][[1]]$Cs, labels2 = mod1$classification, similarity = "jaccard")
  sim.rand.mclust[i] =  cluster_similarity(labels1 = groups[[i]][[1]]$Cs, labels2 = mod1$classification, similarity = "rand")
}

# Figure EDF1M
plot(2:10, sim.rand, type="l", col="dark red", lwd=3, xlab="Number of clusters", ylab="Rand index", ylim=c(0.2, 0.8))
points(2:10, sim.rand.mclust, type="l", col="dark blue", lwd=3,  ylim=c(0.2, 0.8))
