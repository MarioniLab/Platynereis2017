######################################
#### Script to reproduce Figure 3 ####
######################################

setwd("/Users/eling01/Google Drive/")

# Load normalized data
exp.data <- read.table("Platy/6th_approach/Data/Norm/norm_data.txt", sep = "\t")
HVG <- read.table("Platy/6th_approach/Results/HVG.txt", sep = "\t", stringsAsFactors = FALSE)

# Load cluster identities
clusters <- read.table("Platy/6th_approach/Results/stable_clusters_merged.txt", sep = "\t")

# Cell type tree
tree.df <- data.frame(row.names = HVG$GeneNames[HVG$HVG == TRUE], 
                      ANS = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "ANS"]]), 
                      Cilia = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "Cilia"]]),
                      Muscle = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "Muscle"]]),
                      Gut = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "Gut"]]),
                      Trunk = rowMeans(exp.data[HVG$GeneNames[HVG$HVG == TRUE],rownames(clusters)[clusters$clusters == "DiffTrunk"]]))

# Visualize the dendrogram using spearman dissimilarity
dendo <- hclust(as.dist(sqrt((1-cor(log10(tree.df + 1), method = "spearman"))/2)), method = "complete")
plot(dendo)

# Collect and order cells for plotting
cells <- c(rownames(clusters)[clusters$clusters == "ANS"], 
           rownames(clusters)[clusters$clusters == "Gut"],
           rownames(clusters)[clusters$clusters == "Muscle"],
           rownames(clusters)[clusters$clusters == "Cilia"],
           rownames(clusters)[clusters$clusters == "DiffTrunk"])

# Gene name
# Marker genes are: cGMP-PDE_Loc_72981, HNF4, MHC1-4_IB0AAA23AG11EM1_Myo2_KJ405466_KJ834012, "rsph4a_Loc20848", "g587-t1_cirri_antenna_parapodia_palpae_7TM3_GRM7_Contig11280"
gene <- "g587-t1_cirri_antenna_parapodia_palpae_7TM3_GRM7_Contig11280"

# Calculate the Z-score on log10-transformed counts
cells_mean <- mean(as.numeric(log10(exp.data[gene,cells] + 1)))
cells_sd <- sd(log10(as.numeric(exp.data[gene,cells] + 1)))
z <- (log10(exp.data[gene,cells] + 1) - cells_mean)/cells_sd
barplot(as.numeric(z), col = c("red", "yellow", "dark grey", "blue", "green")[clusters[cells,2]])





