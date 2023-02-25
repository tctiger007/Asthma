library(edgeR)
library(DESeq2)
library(readxl)
library(dplyr)
library(stringr)
library(rstudioapi) # only useful in Rstudio
library(gss)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(VennDiagram)
## mute venndiagram's log files 
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

library(gridExtra)



source("./Biofluc_dataLoader.R")
source("./Biofluc_DE.R")
source("./Biofluc_viz.R")


