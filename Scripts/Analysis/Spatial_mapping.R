#########################################################
#### Script for spatial transcriptomics and Figure 2 ####
#########################################################

setwd("/Users/eling01/Google Drive/")

# Helper functions for mapping
source("Platy/6th_approach/Scripts/Analysis/functions_for_plotting.R")

# Load data
mapping <- read.csv("Platy/6th_approach/Data/Spatial_mapping/scaled_results.csv")
coords <- read.csv("Platy/6th_approach/Data/Spatial_mapping/3DCord.csv")
mapping.genes <- read.table("Platy/6th_approach/Data/Spatial_mapping/SupVoxMat_npix7_lowremoved_persig75.txt",sep=" ",header=TRUE)

### Collect data
summary.mapping <- data.frame(row.names = colnames(mapping))

# Collect all mapped coords per cell
mapped.coords <- apply(mapping, 2, function(n){as.data.frame(coords[which(n <= 5),])})

# Calculate volumes of cells
volumes <- unlist(lapply(mapped.coords, nrow))

# Cells that couldn't be mapped
which(volumes == 0)

# Cells that are mapped everywhere
length(which(volumes > 10000))

summary.mapping$volumes <- volumes
plot(log10(summary.mapping$volumes))

# Remove these cells
mapped.coords.1 <- mapped.coords[names(volumes)[volumes < 10000 & volumes > 0]]

library(igraph)

# Find connected voxel clusters per cell
connected.voxels <- lapply(mapped.coords.1, function(i){
  cur_dist <- as.matrix(dist(i, method = "euclidean"))
  cur_dist[lower.tri(cur_dist, diag = TRUE)] <- 100
  inds <- which(cur_dist < 13, arr.ind = TRUE)
  
  graph <- make_graph(edges = apply(inds, 1, function(n){c(rownames(cur_dist)[n[1]], colnames(cur_dist)[n[2]])}))
  graph <- as.undirected(graph)
  clusters(graph)
})

for(i in 1:length(connected.voxels)){
  connected.voxels[[i]]$name <- names(connected.voxels)[i]
}

connected.voxels.1 <- lapply(connected.voxels, function(n){
  cur_clus <- n
  m <- match(cur_clus$membership, which(cur_clus$csize < 4))
  cur_clus$membership <- cur_clus$membership[is.na(m)]
  cur_clus$csize <- table(cur_clus$membership)
  cur_clus$no <- length(unique(cur_clus$membership))
  cur_clus$confidence <- mapping[as.numeric(names(cur_clus$membership)),cur_clus$name]
  cur_clus
})

### Only keep cluster with highest confidence and larger than 5 voxels
connected.voxels.2 <- lapply(connected.voxels.1, function(n){
  cur_clus <- n
  m <- match(cur_clus$membership, which(cur_clus$csize < 5))
  cur_clus$membership <- cur_clus$membership[is.na(m)]
  cur_clus$csize <- table(cur_clus$membership)
  cur_clus$no <- length(unique(cur_clus$membership))
  cur_clus$confidence <- mapping[as.numeric(names(cur_clus$membership)),cur_clus$name]
  cur_clus
})

summary.mapping$clusters <- NA
summary.mapping[names(connected.voxels.1),2] <- unlist(lapply(connected.voxels.1, function(n){n$no}))

# Identify the proportion of each cell mapping to the head or trunk
### Head trunk
head.trunk <- data.frame(row.names = rownames(coords), head = ifelse(coords[,1] < 75, 1, 0), trunk = ifelse(coords[,1] >= 75, 1, 0))

stats <- lapply(connected.voxels.1, function(n){
  cur_stats <- vector(length = 4)
  cur_stats[1] <- length(which(head.trunk[names(n$membership),1] == 1))/length(n$membership)
  cur_stats[2] <- length(which(head.trunk[names(n$membership),2] == 1))/length(n$membership)
  head.count <- 0
  trunk.count <- 0
  for(i in unique(n$membership)){
    if(mean(head.trunk[names(n$membership)[n$membership == i],1]) > 0.5){
      head.count <- head.count + 1
    }
    else{
      trunk.count <- trunk.count + 1
    }
  }
  cur_stats[3] <- head.count
  cur_stats[4] <- trunk.count
  cur_stats
})

summary.mapping$perc.head <- NA
summary.mapping$perc.trunk <- NA
summary.mapping$clus.count.head <- NA
summary.mapping$clus.count.trunk <- NA

summary.mapping[names(stats),3:6] <- matrix(unlist(stats), ncol = 4, byrow = TRUE) 

# Compare mapping to marker genes for each region
# VAT1L, Syt, Pika, PHC2, WNT8, Ank, MHC

ref <- data.frame(row.names = rownames(mapping.genes))
ref$ventral_neuro_Syt <- ifelse(coords[,1] >=75 & mapping.genes[,"Syt"] == 1, 1, 0)
ref$ventral_neuro_Pika <- ifelse(coords[,1] >=75 & mapping.genes[,"Pikachu"] == 1, 1, 0)
ref$apical_neuro_Syt <- ifelse(coords[,1] <75 & mapping.genes[,"Syt"] == 1 & mapping.genes[,"Phc2"] == 0, 1, 0)
ref$apical_neuro_Pika <- ifelse(coords[,1] <75 & mapping.genes[,"Pikachu"] == 1 & mapping.genes[,"Phc2"] == 0, 1, 0)
ref$ANS_Syt <- ifelse(coords[,1] <75 & mapping.genes[,"Syt"] == 1 & mapping.genes[,"Phc2"] == 1, 1, 0)
ref$ANS_Pika <- ifelse(coords[,1] <75 & mapping.genes[,"Pikachu"] == 1 & mapping.genes[,"Phc2"] == 1, 1, 0)
ref$MB_Syt <- ifelse(coords[,1] <75 & mapping.genes[,"Syt"] == 1 & mapping.genes[,"Phc2"] == 0, 1, 0)
ref$MB_Pika <- ifelse(coords[,1] <75 & mapping.genes[,"Pikachu"] == 1 & mapping.genes[,"Phc2"] == 0, 1, 0)
ref$striat_MHC <- ifelse(mapping.genes[,"MHC14"] == 1, 1, 0)
ref$striat_TropI <- ifelse(mapping.genes[,"TropI"] == 1, 1, 0)
ref$cilia <- ifelse(mapping.genes[,"Foxj"] == 1, 1, 0)
ref$NAectoderm <- ifelse(mapping.genes[,"uncx"] == 1, 1, 0)

summary.mapping[colnames(ref)] <- NA

stats.genes <- lapply(connected.voxels.1, function(n){
  cur_stats <- c(0,0,0,0,0,0,0,0,0,0,0,0)
  for(i in unique(n$membership)){
    cur_names <- names(n$membership)[n$membership == i]
    m <- colMeans(ref[cur_names,])
    cur_stats <- cur_stats + ifelse(m > 0.5, 1, 0)
  }
  cur_stats/n$no
})

summary.mapping[names(stats),7:18] <- matrix(unlist(stats.genes), ncol = 12, byrow = TRUE) 

# Compar to cluster identities
clusters <- read.table("Platy/6th_approach/Results/stable_clusters_merged.txt", sep = "\t", header = TRUE)

summary.mapping$prior.clustering <- clusters$clusters[match(rownames(summary.mapping), rownames(clusters))]
summary.mapping <- summary.mapping[order(summary.mapping$prior.clustering),]

summary.mapping <- summary.mapping[!is.na(summary.mapping$perc.head),]

write.csv(summary.mapping, "Platy/6th_approach/Results/mapping_summary.csv")

# Visualize mapping - clusters with highest confidence
# Figure 2 E - Muscle
muscle <- connected.voxels.2[rownames(clusters)[clusters$clusters == "Muscle"]]
muscle <- muscle[!sapply(muscle, is.null)]
muscle <- muscle[sapply(muscle, function(n){n$no}) > 0]
muscle.1 <- lapply(muscle, function(n){
  print(n$name)
  means <- sapply(unique(n$membership), function(x){
    c(mean(coords[names(n$membership)[n$membership == x],1]), 
      mean(coords[names(n$membership)[n$membership == x],2]),
      mean(coords[names(n$membership)[n$membership == x],3]))
  })
  colnames(means) <- unique(n$membership)
  highest <- sapply(unique(n$membership), function(x){
    mean(n$confidence[which(n$membership == x)])
  })
  n$highestConfidence <- names(n$membership[n$membership == as.numeric(colnames(means)[which(highest == min(highest))])])
  n$means <- means[,which(highest == min(highest))]
  n
})

mean.expression.plot(mean.coords = sapply(muscle.1, function(n){n$means}), all.coords = coords, col.var = rep(10, length(muscle.1)),cl = "black", brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"MHC14"] > 0], cluster.col = adjustcolor("dark green", alpha.f = 0.03), theta = 90, phi=90, bty="n")

# All cluster centers
muscle.2 <- lapply(muscle, function(n){
  print(n$name)
  means <- sapply(unique(n$membership), function(x){
    c(mean(coords[names(n$membership)[n$membership == x],1]), 
      mean(coords[names(n$membership)[n$membership == x],2]),
      mean(coords[names(n$membership)[n$membership == x],3]))
  })
  colnames(means) <- unique(n$membership)
  highest <- sapply(unique(n$membership), function(x){
    mean(n$confidence[which(n$membership == x)])
  })
  n$avgConf <- highest
  n$means <- means[,which(highest <= 5)]
  n
})

# Figure 2 F - Example cell
mat <- muscle.2$C11x0501$means
mean.expression.plot(mean.coords = mat, all.coords = coords, col.var = rep(10, ncol(mat)),cl = "blue", brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"MHC14"] > 0], cluster.col = adjustcolor("dark green", alpha.f = 0.03), theta = 90, phi=90, bty="n")

# Find the cluster center with highest confidence for selected cell
mean.expression.plot(mean.coords = sapply(muscle.1, function(n){n$means}), all.coords = coords, col.var = c(9, rep(10, 22)),cl = c("black", "red"), brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"MHC14"] > 0], cluster.col = adjustcolor("dark green", alpha.f = 0.03), theta = 90, phi=90, bty="n")


# Figure 2 A - ANS
ANS <- connected.voxels.2[rownames(clusters)[clusters$clusters == "ANS"]]
ANS <- ANS[!sapply(ANS, is.null)]
ANS <- ANS[sapply(ANS, function(n){n$no}) > 0]
ANS.1 <- lapply(ANS, function(n){
  print(n$name)
  means <- sapply(unique(n$membership), function(x){
    c(mean(coords[names(n$membership)[n$membership == x],1]), 
      mean(coords[names(n$membership)[n$membership == x],2]),
      mean(coords[names(n$membership)[n$membership == x],3]))
  })
  colnames(means) <- unique(n$membership)
  highest <- sapply(unique(n$membership), function(x){
    mean(n$confidence[which(n$membership == x)])
  })
  n$highestConfidence <- names(n$membership[n$membership == as.numeric(colnames(means)[which(highest == min(highest))])])
  n$means <- means[,which(highest == min(highest))[1]]
  n
})

mean.expression.plot(mean.coords = sapply(ANS.1, function(n){n$means}), all.coords = coords, col.var = rep(10, length(ANS.1)),cl = "black", brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"Phc2"] > 0], cluster.col = adjustcolor("red", alpha.f = 0.03), theta = 90, phi=90, xlab="Medio-lateral", ylab="Posterior-Anterior", zlab="Dorso-Ventral", bty="n")

# All voxel centres
ANS.2 <- lapply(ANS, function(n){
  print(n$name)
  means <- sapply(unique(n$membership), function(x){
    c(mean(coords[names(n$membership)[n$membership == x],1]), 
      mean(coords[names(n$membership)[n$membership == x],2]),
      mean(coords[names(n$membership)[n$membership == x],3]))
  })
  colnames(means) <- unique(n$membership)
  highest <- sapply(unique(n$membership), function(x){
    mean(n$confidence[which(n$membership == x)])
  })
  n$avgConf <- highest
  n$means <- means[,which(highest <= 5)]
  n
})

# Figure 2 B - Example cell
mat <- ANS.2$C31x8101$means
mean.expression.plot(mean.coords = mat, all.coords = coords, col.var = rep(10, ncol(mat)),cl = "blue", brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"Phc2"] > 0], cluster.col = adjustcolor("red", alpha.f = 0.03), theta = 90, phi=90, bty="n")

# Find the cluster center with highest confidence for selected cell
mean.expression.plot(mean.coords = sapply(ANS.1, function(n){n$means}), all.coords = coords, col.var = c(rep(10, 4), 9, rep(10, 18)),cl = c("black", "red"), brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"Phc2"] > 0], cluster.col = adjustcolor("red", alpha.f = 0.03), theta = 90, phi=90, xlab="Medio-lateral", ylab="Posterior-Anterior", zlab="Dorso-Ventral", bty="n")

# Figure 2 I - Non apical ectoderm
Trunk <- connected.voxels.2[rownames(clusters)[clusters$clusters == "DiffTrunk"]]
Trunk <- Trunk[!sapply(Trunk, is.null)]
Trunk <- Trunk[sapply(Trunk, function(n){n$no}) > 0]
Trunk.1 <- lapply(Trunk, function(n){
  print(n$name)
  means <- sapply(unique(n$membership), function(x){
    c(mean(coords[names(n$membership)[n$membership == x],1]), 
      mean(coords[names(n$membership)[n$membership == x],2]),
      mean(coords[names(n$membership)[n$membership == x],3]))
  })
  colnames(means) <- unique(n$membership)
  highest <- sapply(unique(n$membership), function(x){
    mean(n$confidence[which(n$membership == x)])
  })
  n$highestConfidence <- names(n$membership[n$membership == as.numeric(colnames(means)[which(highest == min(highest))])])
  n$means <- means[,which(highest == min(highest))[1]]
  n
})

mean.expression.plot(mean.coords = sapply(Trunk.1, function(n){n$means}), all.coords = coords, col.var = rep(10, length(Trunk.1)),cl = "black", brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"uncx"] > 0], cluster.col = adjustcolor("black", alpha.f = 0.03), theta = 90, phi=90, xlab="Medio-lateral", ylab="Posterior-Anterior", zlab="Dorso-Ventral", bty="n")

# All voxel clusters
Trunk.2 <- lapply(Trunk, function(n){
  print(n$name)
  means <- sapply(unique(n$membership), function(x){
    c(mean(coords[names(n$membership)[n$membership == x],1]), 
      mean(coords[names(n$membership)[n$membership == x],2]),
      mean(coords[names(n$membership)[n$membership == x],3]))
  })
  colnames(means) <- unique(n$membership)
  highest <- sapply(unique(n$membership), function(x){
    mean(n$confidence[which(n$membership == x)])
  })
  n$avgConf <- highest
  n$means <- means[,which(highest <= 5)]
  n
})

# Figure 2 J - Plot example cell
mat <- Trunk.2$C11x0901$means
mean.expression.plot(mean.coords = mat, all.coords = coords, col.var = rep(10, ncol(mat)),cl = "blue", brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"uncx"] > 0], cluster.col = adjustcolor("black", alpha.f = 0.03), theta = 90, phi=90, bty="n")

# Find the cluster center with highest confidence for selected cell
mean.expression.plot(mean.coords = sapply(Trunk.1, function(n){n$means}), all.coords = coords, col.var = c(9, rep(10, length(Trunk.1)-1)),cl = c("black", "red"), brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"uncx"] > 0], cluster.col = adjustcolor("black", alpha.f = 0.03), theta = 90, phi=90, xlab="Medio-lateral", ylab="Posterior-Anterior", zlab="Dorso-Ventral", bty="n")

# Figure 2 G - Cilia
Cilia <- connected.voxels.2[rownames(clusters)[clusters$clusters == "Cilia"]]
Cilia <- Cilia[!sapply(Cilia, is.null)]
Cilia.1 <- lapply(Cilia, function(n){
  print(n$name)
  means <- sapply(unique(n$membership), function(x){
    c(mean(coords[names(n$membership)[n$membership == x],1]), 
      mean(coords[names(n$membership)[n$membership == x],2]),
      mean(coords[names(n$membership)[n$membership == x],3]))
  })
  colnames(means) <- unique(n$membership)
  highest <- sapply(unique(n$membership), function(x){
    mean(n$confidence[which(n$membership == x)])
  })
  n$highestConfidence <- names(n$membership[n$membership == as.numeric(colnames(means)[which(highest == min(highest))])])
  n$means <- means[,which(highest == min(highest))[1]]
  n
})

mean.expression.plot(mean.coords = sapply(Cilia.1, function(n){n$means}), all.coords = coords, col.var = rep(10, length(Cilia.1)),cl = "black", brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"Foxj"] > 0], cluster.col = adjustcolor("yellow", alpha.f = 0.1), theta = 90, phi=90, xlab="Medio-lateral", ylab="Posterior-Anterior", zlab="Dorso-Ventral", bty="n")

### Threshold on cluster centres
Cilia.2 <- lapply(Cilia, function(n){
  print(n$name)
  means <- sapply(unique(n$membership), function(x){
    c(mean(coords[names(n$membership)[n$membership == x],1]), 
      mean(coords[names(n$membership)[n$membership == x],2]),
      mean(coords[names(n$membership)[n$membership == x],3]))
  })
  colnames(means) <- unique(n$membership)
  highest <- sapply(unique(n$membership), function(x){
    mean(n$confidence[which(n$membership == x)])
  })
  n$avgConf <- highest
  n$means <- means[,which(highest <= 5)]
  n
})

# Figure 2 H - Plot example cells
mat <- Cilia.2$C22x0201L$means
mean.expression.plot(mean.coords = mat, all.coords = coords, col.var = rep(10, ncol(mat)),cl = "blue", brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"Foxj"] > 0], cluster.col = adjustcolor("yellow", alpha.f = 0.03), theta = 90, phi=90, bty="n")

# Find the cluster center with highest confidence for selected cell
mean.expression.plot(mean.coords = sapply(Cilia.1, function(n){n$means}), all.coords = coords, col.var = c(9, rep(10, length(Cilia.1) - 1)),cl = c("black", "red"), brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"Foxj"] > 0], cluster.col = adjustcolor("yellow", alpha.f = 0.1), theta = 90, phi=90, xlab="Medio-lateral", ylab="Posterior-Anterior", zlab="Dorso-Ventral", bty="n")

# Figure 2 C - Gut cells
### Plot Gut cells
Gut <- connected.voxels.2[rownames(clusters)[clusters$clusters == "Gut"]]
Gut <- Gut[!sapply(Gut, is.null)]
Gut.1 <- lapply(Gut, function(n){
  print(n$name)
  means <- sapply(unique(n$membership), function(x){
    c(mean(coords[names(n$membership)[n$membership == x],1]), 
      mean(coords[names(n$membership)[n$membership == x],2]),
      mean(coords[names(n$membership)[n$membership == x],3]))
  })
  colnames(means) <- unique(n$membership)
  highest <- sapply(unique(n$membership), function(x){
    mean(n$confidence[which(n$membership == x)])
  })
  n$highestConfidence <- names(n$membership[n$membership == as.numeric(colnames(means)[which(highest == min(highest))])])
  n$means <- means[,which(highest == min(highest))[1]]
  n
})

mean.expression.plot(mean.coords = sapply(Gut.1, function(n){n$means}), all.coords = coords, col.var = rep(10, length(Gut.1)),cl = "black", brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"HNF4"] > 0], cluster.col = adjustcolor("blue", alpha.f = 0.1), theta = 90, phi=90, xlab="Medio-lateral", ylab="Posterior-Anterior", zlab="Dorso-Ventral", bty="n")

### Threshold on cluster centres
Gut.2 <- lapply(Gut, function(n){
  print(n$name)
  means <- sapply(unique(n$membership), function(x){
    c(mean(coords[names(n$membership)[n$membership == x],1]), 
      mean(coords[names(n$membership)[n$membership == x],2]),
      mean(coords[names(n$membership)[n$membership == x],3]))
  })
  colnames(means) <- unique(n$membership)
  highest <- sapply(unique(n$membership), function(x){
    mean(n$confidence[which(n$membership == x)])
  })
  n$avgConf <- highest
  n$means <- means[,which(highest <= 5)]
  n
})

# Figure 2 D - Plot example cell
mat <- Gut.2$C5x2301L$means
mean.expression.plot(mean.coords = mat, all.coords = coords, col.var = rep(10, ncol(mat)),cl = "blue", brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"HNF4"] > 0], cluster.col = adjustcolor("blue", alpha.f = 0.03), theta = 90, phi=90, bty="n")

# Find the cluster center with highest confidence for selected cell
mean.expression.plot(mean.coords = sapply(Gut.1, function(n){n$means}), all.coords = coords, col.var = c(10,9,10,10),cl = c("black", "red"), brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"HNF4"] > 0], cluster.col = adjustcolor("blue", alpha.f = 0.1), theta = 90, phi=90, xlab="Medio-lateral", ylab="Posterior-Anterior", zlab="Dorso-Ventral", bty="n")

#### Extended Data Figure 5
# Find cells mapping to Syt1 region in head
Syt1.head <- lapply(ANS.1, function(n){
  n$Syt <- ifelse(sum(mapping.genes[n$highestConfidence,"Syt"])/length(n$highestConfidence) > 0.5, TRUE, FALSE)
  n
})

mean.expression.plot(mean.coords = sapply(ANS.1, function(n){n$means}), all.coords = coords, col.var = as.numeric(unlist(lapply(Syt1.head, function(n){n$Syt})))+1,cl = c("black", "red"), brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"Syt"] > 0], cluster.col = adjustcolor("purple", alpha.f = 0.03), theta = 90, phi=90, xlab="Medio-lateral", ylab="Posterior-Anterior", zlab="Dorso-Ventral", bty="n")

# Find cells mapping to Syt1 in non-apical surface
Syt1.trunk <- lapply(Trunk.1, function(n){
  n$Syt <- ifelse(sum(mapping.genes[n$highestConfidence,"Syt"])/length(n$highestConfidence) > 0.5, TRUE, FALSE)
  n
})

mean.expression.plot(mean.coords = sapply(Trunk.1, function(n){n$means}), all.coords = coords, col.var = as.numeric(unlist(lapply(Syt1.trunk, function(n){n$Syt})))+1,cl = c("black", "red"), brain.cl = adjustcolor("grey", alpha.f = 0.01), col.var.title = "test", add.cluster = rownames(coords)[mapping.genes[,"Syt"] > 0], cluster.col = adjustcolor("purple", alpha.f = 0.03), theta = 90, phi=90, xlab="Medio-lateral", ylab="Posterior-Anterior", zlab="Dorso-Ventral", bty="n")
