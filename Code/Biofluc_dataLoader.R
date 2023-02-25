library(rstudioapi)
library(DESeq2)
library(readxl)
library(dplyr)
library(stringr)
library(gss)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(VennDiagram)
## mute venndiagram's log files 
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
library(gridExtra)


path = dirname(getActiveDocumentContext()$path)
setwd(path)

Counts = read.csv("../Data/Run268to294.csv", check.names = FALSE)
colnames(Counts)[1] = "miRNAs"
rownames(Counts) = Counts$miRNAs
Counts = Counts[, -1]
Counts = Counts + 1e-8


################################################################################
######################################## PCA ###################################
################################################################################
Traits = read.csv("../Data/Run268to294traits.csv", check.names = FALSE)  
colnames(Traits)[which(colnames(Traits) == "Days_from_Exposure")] = "Time"


Traits$Asthma = ifelse(Traits$Subject_Group == "H", 0, 1)
Traits[, c("Run", "Subject_ID", "Asthma")] = 
  lapply(Traits[, c("Run", "Subject_ID", "Asthma")], factor)

CountsWithGroupInfo = cbind.data.frame(Traits, t(Counts))
rownames(CountsWithGroupInfo) = CountsWithGroupInfo$Sequence_ID

################################################################################
################# preprocess: filtering genes and samples ######################
######################## preprocess: Normalization #############################
################################################################################
### remove outliers (samples) based on the lib size (filter samples)
lib.size = data.frame(colSums(Counts))
quantile(lib.size$colSums.Counts., seq(0, 1, by = .01))

remove.samples = rownames(lib.size)[which(lib.size$colSums.Counts.< 4232
                                          |lib.size$colSums.Counts. > 188707)]  
remove.samples = unique(remove.samples) 
colSums(Counts)
df = cbind.data.frame(Traits, colSums(Counts))

Counts = Counts[,-which(colnames(Counts) %in% c(remove.samples))]
CountsWithGroupInfo = CountsWithGroupInfo[-which(CountsWithGroupInfo$Sequence_ID %in% c(remove.samples)),]
Traits = Traits[-which(Traits$Sequence_ID %in% c(remove.samples)),]

###### filter out genes that have low counts 
keep = rowSums(Counts)> 200
Counts = Counts[keep,]
Counts = as.matrix(Counts)
CountsWithGroupInfo = cbind.data.frame(CountsWithGroupInfo[,c(1:11)],
                                       t(Counts))

### Normalization deseq2
cData <- data.frame(Subject_Group = as.factor(CountsWithGroupInfo$Subject_Group),
                    Run = as.factor(CountsWithGroupInfo$Run),
                    Time = CountsWithGroupInfo$Time,
                    InfectionStatus = CountsWithGroupInfo$InfectionStatus)
rownames(cData) = CountsWithGroupInfo$Sequence_ID
Counts.r = round(Counts)
dds = DESeqDataSetFromMatrix(countData = Counts.r, colData = cData, design = ~ 1)
dds = estimateSizeFactors(dds)
# DESeq2 normalization counts
Counts.normalized.deseq2 = counts(dds, normalized = TRUE)
CountsWithGroupInfo.normalized.deseq2 = cbind.data.frame(CountsWithGroupInfo[,c(1:11)],
                                                         t(Counts.normalized.deseq2))
# write.csv(CountsWithGroupInfo.normalized, "CountsWithGroupInfo_normalized_forgss.csv")
CountsWithGroupInfo.normalized = CountsWithGroupInfo.normalized.deseq2
# write.csv(CountsWithGroupInfo.normalized, "CountsWithGroupInfo_normalized_forgss.csv")
# saveRDS(CountsWithGroupInfo.normalized, "CountsWithGroupInfo_normalized.rds")

CountsWithGroupInfo.normalized = CountsWithGroupInfo.normalized[,-1]
factor.idx = which(colnames(CountsWithGroupInfo.normalized) %in% c("Subject_ID","Subject_Group","Visit","TimePoint", 
                                                                   "Run", "InfectionStatus", "SampleBins", "Asthma"))
CountsWithGroupInfo.normalized = cbind.data.frame(CountsWithGroupInfo.normalized[,-factor.idx], 
                                                  lapply(CountsWithGroupInfo.normalized[, factor.idx],
                                                         as.factor))
CountsWithGroupInfo.normalized = CountsWithGroupInfo.normalized[,c(1:3, 
                                                                   (dim(CountsWithGroupInfo.normalized)[2]-6):dim(CountsWithGroupInfo.normalized)[2],
                                                                   4: (dim(CountsWithGroupInfo.normalized)[2]-7))]


#### clinical data including spirometry (first measurement) data 
traitData1 = read.csv("../Data/MactchedClinicalDataMeasure1.csv",
                      check.names = FALSE, row.names = 1)
traitData1$Sequence_ID = paste0(str_remove(traitData1$Subject_ID, "^0+"), "V",
                                traitData1$Visit)
dim(traitData1) #241  50
names(traitData1) 
colnames(traitData1)[which(colnames(traitData1) == "FEV1/FVC")] = "FEV1FVC"

commonID = traitData1$Sequence_ID[traitData1$Sequence_ID %in% rownames(CountsWithGroupInfo.normalized)]
# keep common IDs (in both traits and sequencing data)
traitData2 = traitData1[match(commonID,traitData1$Sequence_ID),]
CountsWithGroupInfo.normalized1 = CountsWithGroupInfo.normalized[match(commonID,rownames(CountsWithGroupInfo.normalized)),]

df = cbind.data.frame(traitData2, CountsWithGroupInfo.normalized1)
df = df[, -c(1:9, 49:60)]

rownames(df) = traitData2$Sequence_ID

removetraits = which(colnames(df) %in% c("NL_vol", "Weight_gel_NL", "DTT", "DNAse", "cells_total",
                                         "cell_and_epi", "dilution", "cells"))
df = df[, -removetraits]

Traits1 = Traits[match(commonID, Traits$Sequence_ID),]

CountsWithGroupInfo.normalized1 = CountsWithGroupInfo.normalized[, -c(1:10)]
CountsWithGroupInfo.normalized1 = CountsWithGroupInfo.normalized1[match(commonID, rownames(CountsWithGroupInfo.normalized1)),]

Data = list(CountsWithGroupInfo.normalized1, 
            Traits1,
            traitData2)

names(Data) = c("Counts", "Meta", "ClinicalTraits")

# save(Data, file = "./data.RData")


