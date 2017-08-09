######################################################
#### Script to estimate technical noise per batch ####
######################################################

setwd("/Users/eling01/Google Drive/")

library(BASiCS)

# Read in data
raw.reads <- read.table("Platy/6th_approach/Data/Raw/QC_rawData.tsv", sep = "\t")

# Read in ERCC concentrations
ERCC.conc <- read.table("Platy/6th_approach/Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

# Calculate number of spike-ins
ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))
ERCC.num.final <- ERCC.num/800000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

# Separate cells by batches
chips <- sapply(colnames(raw.reads), function(n){paste(unlist(strsplit(n, split = "x"))[1], "x", sep = "")})

# Lists to store BASiCS Data and MCMC objects
MCMCs <- list()
Datas <- list()

for(i in 2:length(unique(chips))){
  # Select cells of current batch
  Counts <- raw.reads[,grepl(unique(chips)[i], colnames(raw.reads))]
  Counts <- Counts[rowMeans(Counts) > 1,]
  
  # Prepare BASiCS Data object
  Tech <- grepl("ERCC-0", rownames(Counts))
  SpikeInput <- ERCC.num.final[rownames(Counts)[grepl("ERCC-0", rownames(Counts))],1]
  SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                             "Molecules" = SpikeInput,
                             stringsAsFactors = FALSE)
  
  Data <- newBASiCS_Data(as.matrix(Counts), Tech, SpikeInput.1)
  
  MCMC_Output <- BASiCS_MCMC(Data, N = 20000, Thin = 20, Burn = 10000, PrintProgress = TRUE, PriorDelta = 'log-normal')
  
  MCMCs[[paste("Batch", unique(chips)[i])]] <- MCMC_Output
  Datas[[paste("Batch", unique(chips)[i])]] <- Data
}

rm(list=setdiff(ls(), c("MCMCs", "Datas")))

save.image("Platy/6th_approach/Results/batches_MCMC.RData")

BASiCS_VarianceDecomp(Data = Datas$`Batch C1x`, object = MCMCs$`Batch C1x`)
BASiCS_VarianceDecomp(Data = Datas$`Batch C5x`, object = MCMCs$`Batch C5x`)
BASiCS_VarianceDecomp(Data = Datas$`Batch C4x`, object = MCMCs$`Batch C4x`)
BASiCS_VarianceDecomp(Data = Datas$`Batch C31x`, object = MCMCs$`Batch C31x`)
BASiCS_VarianceDecomp(Data = Datas$`Batch C26x`, object = MCMCs$`Batch C26x`)
BASiCS_VarianceDecomp(Data = Datas$`Batch C25x`, object = MCMCs$`Batch C25x`)
BASiCS_VarianceDecomp(Data = Datas$`Batch C22x`, object = MCMCs$`Batch C22x`)
BASiCS_VarianceDecomp(Data = Datas$`Batch C11x`, object = MCMCs$`Batch C11x`)
BASiCS_VarianceDecomp(Data = Datas$`Batch CN61x`, object = MCMCs$`Batch CN61x`)
