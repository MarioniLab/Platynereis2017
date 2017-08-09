######################################################
#### Script to generate large heatmap of figure 4 ####
######################################################

setwd("/Users/eling01/Google Drive/")

library(pheatmap)

# Load normalized data
exp.data <- read.table("Platy/6th_approach/Data/Norm/norm_data.txt", sep = "\t")

# Read in highly variable genes
HVG <- read.table("Platy/6th_approach/Results/HVG.txt", sep = "\t", stringsAsFactors = FALSE)

# Load cluster identities
clusters <- read.table("Platy/6th_approach/Results/stable_clusters.txt", sep = "\t")

# Load gene annotation
gene.names <- read.table(file = "/Users/eling01/Dropbox/Platy_neurons/Results/allgenes_abbr.txt", sep = "\t", header = TRUE, fill = TRUE, quote = "", stringsAsFactors = FALSE)

# Selected group specifc genes for displaying
ANS.genes <- c("V9I5J3_TIMP1_Loc_5954", "PHC2", "GNXQN_neuropeptide_precursor", "Syt17", "Syt4", "WI_neuropeptide_precursor", "Pdu_5HTR1A", "GGNamide", "Loc_6730_VAT1L", "Contig1452_Neuroglobin", "MIP_allatostatin_B_neuropeptide_precursor", "allatostatin-C", "insulin_related_peptide_2_precursor",   "Insulin-like_D1-2_g30672", "NPY_Primr", "vasotocin-neurophysin", "Peropsin",  "cGMP-PDE_Loc_72981", "CNGAalpha", "CNGB", "Pigment_dispersing_factor", "proenkephalin", "PDE-5A_438", "cOpsin",  "c-opsin2", "NK5", "FoxQ2", "Bsx","TBR1_MOUSE T-box brain protein 1")
Gut.genes <- c("Sytalpha", "PLCdelta-like_405", "FMRF-R_Loc_27555", "Contactin_Loc_10492", "CCWamide_neuropeptide_precursor", "ChymotrypsinA_Loc_91561_tryptase-gamma", "PY_neuropeptide_precursor", "WLDpeptide", "HNF4")
Muscle.genes <- c("FilaminA_1K_75_Loc_1862Ischia", "mrlc_Contig17840", "Pdu-troponinT", "MHC1-4_IB0AAA23AG11EM1_Myo2_KJ405466_KJ834012","P70566_tropmodulin_1A_75_Loc_1964Ischia", "Pdu_Col4a2_Contig4086", "Smoothelin", "RyanodinR2_6K_82_Loc_365Ischia", "RYR44F_RYR2_ENR2", "MyoD", "HAND",  "Mox", "Pitx",   "Paraxis_scleraxis")
Trunk.genes <- c("Esyt2a", "MMEP_Loc_112633", "LAMIN_N2-15_g103177", "NF70_LOLPE 70 kDa neurofilament protein", "Pdu_trscr_assembly_1|91867", "Msx", "ETS4_DROME DNA-binding protein D-ETS-4", "Sema2a_R20", "R7UV80_Anoctamin_Loc_432", "A0A067R0I3_PCDH15_Loc_34801", "g587-t1_cirri_antenna_parapodia_palpae_7TM3_GRM7_Contig11280", "CD109_Contig10163", "Q6ITT8_CathepsinL_Loc_47202", "grainyhead_comp408863_Loc_75846", "EHF_BOVIN ETS homologous factor", "nachr-alpha6_7_9_455", "nat2", "cand.HSPG2_Ischia", "Pdu_trscr_assembly_1|102268", "NRT_DROME Neurotactin", "Pdu_similar-to-neurotrypsin", "NETR_PONPY Neurotrypsin", "Pdu_trscr_assembly_1|3542", "Pdu_trscr_assembly_1|15925", "SLC22", "VEGFR_FGFR_C1-13", "g47357-t1_palpae_specific_g52479-t1_antenna_cirri_palpae_iGlur_Contig10681", "Pdu_Slc1a2_putative_alt_5prim-UTR", "Hes1","Hes2", "Hes6",  "Hes11", "Wnt4", "uncx4", "Phox2", "Tbx20_Loc_14134", "Hb9_Mnx")
Cilia.genes <- c("rsph1_Loc36539",  "rsph4a_Loc20848", "Tektin-4", "K1R508_IFT80_6K_82_Loc_1122Ischia",  "Q62559_IFT52_Loc_27660", "K1Q768_IFT57_6K_82_Loc_15712", "Q6NWV3_IFT122_Contig3098",   "Q61025_IFT20_Loc_42030", "K1QZL1_IFT88_6A_75_Loc_9776Ischia", "K1Q5N1_Spag8_Loc_73606", "R7TH40_Kinesin-like_Loc_31446","foxJ", "Gbx",  "Tektin-3",  "tailless", "RX") 

# Find gene names for annotation
annot.names <- gene.names[match(c(ANS.genes, Gut.genes,Muscle.genes,Trunk.genes, Cilia.genes), gene.names[,1]),2]

# Match with Jekely annotation
GJ.genes <- read.table("Platy/6th_approach/Data/Gene_annotation/names_GJ_split.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
GJ.genes$Evalue <- as.numeric(sapply(sapply(as.character(GJ.genes$cut_blasthit), function(n){unlist(strsplit(n, split = "E-value: "))[2]}), function(x){unlist(strsplit(x, split = " "))[1]}))
GJ.genes$BLAST_hit <- as.character(sapply(sapply(as.character(GJ.genes$cut_blasthit), function(n){unlist(strsplit(n, split = "\\|"))[3]}), function(x){unlist(strsplit(x, split = " OS"))[1]}))

annot.names[is.na(annot.names)] <- gene.names[match(GJ.genes$jekely_ID[match(c(ANS.genes, Gut.genes,Muscle.genes,Trunk.genes, Cilia.genes)[which(is.na(annot.names))], GJ.genes$BLAST_hit)], gene.names[,1]),2]

# Plot heatmap
pheatmap(log10(exp.data[c(ANS.genes, Gut.genes,Muscle.genes,Trunk.genes, Cilia.genes),rownames(clusters)[-which(grepl("Undiff", clusters$clusters))]] + 1), 
         cluster_rows = FALSE, 
         annotation_col = data.frame(row.names = rownames(clusters)[-which(grepl("Undiff", clusters$clusters))], group = clusters$clusters[-which(grepl("Undiff", clusters$clusters))]), 
         clustering_distance_cols = as.dist((1- cor(log10(exp.data[as.character(HVG$GeneNames[HVG$HVG == TRUE]), rownames(clusters)[-which(grepl("Undiff", clusters$clusters))]] + 1), method = "spearman"))/2), 
         col = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100), 
         cellwidth = 4, 
         cellheight = 4, 
         gaps_row = c(29, 38, 52, 89), 
         fontsize = 4,
         labels_row = annot.names)

