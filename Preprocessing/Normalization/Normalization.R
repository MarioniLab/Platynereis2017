#####################################################################
#### Script to normalize the quality filtered raw counts data ######
#####################################################################

library(BASiCS)

# Read in the raw data
raw.data <- read.table("/Users/eling01/Google Drive/Platy/6th_approach/Data/Raw/QC_rawData.tsv", sep = "\t")

# Read in the ERCC concentrations
ERCC.conc <- read.table("/Users/eling01/Google Drive/Platy/6th_approach/Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

# Calculate transcript count per spike-in 
ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))
ERCC.num.final <- ERCC.num/800000 # Dilution factor
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

# Prepare BASiCS object
Tech <- grepl("ERCC-0", rownames(raw.data))
SpikeInput <- ERCC.num.final[rownames(raw.data)[grepl("ERCC-0", rownames(raw.data))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput,
                           stringsAsFactors = FALSE)

# Include batch information
chips <- sapply(colnames(raw.data), function(n){unlist(strsplit(n, split = "x"))[1]})

Data <- newBASiCS_Data(as.matrix(raw.data), Tech, SpikeInput.1, BatchInfo = chips)

# Run the MCMC
MCMC_Output <- BASiCS_MCMC(Data, N = 40000, Thin = 20, Burn = 20000, PrintProgress = TRUE)

# Check for convergence of chains
library(coda)
plot(mcmc(MCMC_Output@mu[,1]))
plot(mcmc(MCMC_Output@delta[,1]))
plot(mcmc(MCMC_Output@phi[,1]))
plot(mcmc(MCMC_Output@s[,1]))
plot(mcmc(MCMC_Output@nu[,1]))
plot(mcmc(MCMC_Output@theta[,1]))

# Calculate denoised counts
DenoisedCounts = BASiCS_DenoisedCounts(Data = Data, Chain = MCMC_Output)
exp.data <- DenoisedCounts

# Calculate highly variable genes
HVG <- BASiCS_DetectHVG(Data = Data, object = MCMC_Output, VarThreshold = 0.98, EviThreshold = 0.7, OrderVariable = "Mu")

rm(list=setdiff(ls(), c("exp.data", "HVG", "MCMC_Output", "Data")))

save.image("/Users/eling01/Google Drive/Platy/6th_approach/Data/Norm/norm_data.RData")
write.table(as.data.frame(exp.data), "/Users/eling01/Google Drive/Platy/6th_approach/Data/Norm/norm_data.txt", sep = "\t")
