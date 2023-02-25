library(Mfuzz)
library(rstudioapi)
library(stringr)
library(Biobase)
library(RColorBrewer)
library(gplots)
library(dplyr)
library(readxl)
library(gss)
library(reshape2)
library(VennDiagram)
library(gridExtra)
library(grid)
library(xlsx)
library(plot3D)
library(fdapace)
library(ggplot2)
library(png)
library(fields)

set.seed(1)
path = dirname(getActiveDocumentContext()$path)
setwd(path)
fit.df = read.csv("../Data/fitted_DEmiRNA.csv")
DE.group3 = read_excel("../Data/sig_gss.xlsx", sheet = 3)
DE.group3 = DE.group3$miRNAs

# # --------------------------------- 2. use lfc ---------------------------------
# Group 3 is the miRNAs that are differentially expressed between asthmatic and healthy subjects at post-RV
source("fit_expr.R")
expr.DE.m3 = fit.expr(DE.group3, log = FALSE)[[2]]
lfc.3 = log2(expr.DE.m3[,c(1:20)]/expr.DE.m3[,c(21:40)])
colnames(lfc.3) = paste0("lfc", seq(1,20,1))
lfc.3.post = lfc.3[, c(15:20)]

color = colorRampPalette(brewer.pal(8, "RdYlGn"))(25)
heatmap.2(lfc.3.post, scale = "column", col = color, key = TRUE,
          trace = "none", density.info = "none",
          key.title = NA,
          lmat=rbind(c(0,3),c(2,1),c(0,4)),
          lwid = c(1.5,4), lhei = c(1.5,4,1.2))


pData = colnames(lfc.3.post)
pData = as.data.frame(pData)
pData$Bin = gsub("[aA-zZ]+", "", pData$pData)
rownames(pData) = pData$pData
pData = subset(pData, select = "Bin")

all(rownames(pData)==colnames(lfc.3.post))
metadata = data.frame(labelDescription="Post-RV Bins",
                      row.names=colnames(pData))
phenoData = new("AnnotatedDataFrame",
                 data=pData, varMetadata=metadata)
head(pData(phenoData))
annotation = rownames(lfc.3.post)
experimentData = new("MIAME",
                      name="WW",
                      lab="FP Lab",
                      title="miRNA")

lfc.3.post.m = as.matrix(lfc.3.post)
exampleSet <- ExpressionSet(assayData=lfc.3.post.m,
                            phenoData=phenoData,
                            experimentData=experimentData,
                            annotation=rownames(lfc.3.post))
exampleSet.r = filter.NA(exampleSet)
exampleSet.s = standardise(exampleSet.r)
mestimate(exampleSet.s)   #3.589289


cselection(exampleSet.s,3.589289,crange=seq(2,18,1),repeats=5,visu=TRUE)

Dmin(exampleSet.s, m=3.589289, crange=seq(2,10,1), repeats=3, visu=TRUE)

### edit function mfuzz.plot to customize the plot 
# trace(mfuzz.plot, edit = TRUE)
# time.labels = as.character(round(seq(-71,27, length.out = 20))),
colo <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700", 
                   "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00", 
                   "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", 
                   "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", 
                   "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7", 
                   "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF", 
                   "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", 
                   "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF", 
                   "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF", 
                   "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", 
                   "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", 
                   "#FF0060", "#FF0048", "#FF0030", "#FF0018")
                   #### 2 clusters
cl = mfuzz(exampleSet.s, c = 2, m = 3.589289)
layout(t(1:3), widths = c(6,6,0.5))
par(mar=rep(2, 4), oma=rep(1, 4), las=1)

# png("./mfuzz/asthma_vs_healthy_post_2clusters_usinglfc.png", width = 12, height = 5, units = "in", res = 300)
mfuzz.plot(exampleSet.s, cl = cl, mfrow = c(1,2),
           time.labels = round(seq(-71,27, length.out = 20)[15:20]),
           new.window = FALSE)
image(1, seq(0, 256, length.out = length(colo)), t(seq_along(colo)), col=colo, axes=FALSE)
# ggarrange()

# dev.off()

membership = cl[["membership"]]
class(membership)
# View(membership)
colnames(membership) = c("Cluster1", "Cluster2") #, "Cluster3")

color = colorRampPalette(brewer.pal(8, "PiYG"))(25)
# check how similar the clusters are
# png("./mfuzz/asthma_vs_healthy_mfuzz_2clusters_usinglfc_membership.png", width = 6, height = 6, units = "in", res = 300)
heatmap.2(membership, scale = "column", col = color, key = TRUE,
          trace = "none", density.info = "none",
          key.title = NA, margins = c(8, 8),
          lmat=rbind(c(0,3),c(2,1),c(0,4)),
          lwid = c(1,4), lhei = c(1.5,4,1.2))
# dev.off()

# use membership > 0.55 as threshold for defining sig miRNAs
summary(membership)
# 3rd Qu.:0.6672   3rd Qu.:0.6809  
# Max.   :0.7766   Max.   :0.7609  
membership = as.data.frame(membership)
c1 = rownames(membership)[which(membership$Cluster1 > 0.55)] #12
c2 = rownames(membership)[which(membership$Cluster2 > 0.55)] #11


## mean function: PACE
lfc = exprs(exampleSet.s)
lfc.melt = melt(lfc)
n = dim(membership)[1]
lfc.melt$Time = c(rep(1,n), rep(6,n), rep(12,n), rep(17,n), rep(22,n), rep(27,n))
lfc.melt$miRNA = c(rep(seq(1,n,1),6))
colnames(lfc.melt)[3] = "logFC(asthmatic/healthy)"

lfc.melt.c1 = lfc.melt %>% filter(Var1 %in% c1)
lfc.melt.c2 = lfc.melt %>% filter(Var1 %in% c2)

Files.c1 =  MakeFPCAInputs(lfc.melt.c1$miRNA, lfc.melt.c1$Time, lfc.melt.c1$`logFC(asthmatic/healthy)`)
Files.c2 =  MakeFPCAInputs(lfc.melt.c2$miRNA, lfc.melt.c2$Time, lfc.melt.c2$`logFC(asthmatic/healthy)`)

fpcaObjFiles.c1 <- FPCA(Files.c1$Ly, Files.c1$Lt, list(plot = TRUE, methodMuCovEst = 'smooth', userBwCov = 2))
fpcaObjFiles.c2 <- FPCA(Files.c2$Ly, Files.c2$Lt, list(plot = TRUE, methodMuCovEst = 'smooth', userBwCov = 2))

fpca.x.c1 = fpcaObjFiles.c1[["workGrid"]]; fpca.y.c1 = fpcaObjFiles.c1[["mu"]]
fpca.x.c2 = fpcaObjFiles.c2[["workGrid"]]; fpca.y.c2 = fpcaObjFiles.c2[["mu"]]

# fpca x and y for plotting 
fpca.c1 = data.frame(x = fpca.x.c1, y = fpca.y.c1)
fpca.c2 = data.frame(x = fpca.x.c2, y = fpca.y.c2)

membership.c1 = membership[match(c1, rownames(membership)),]
membership.c2 = membership[match(c2, rownames(membership)),]

ordered.c1 = rownames(membership.c1)[sort.list(membership.c1$Cluster1, decreasing = T)]
ordered.c2 = rownames(membership.c2)[sort.list(membership.c2$Cluster2, decreasing = T)]

lfc.c1.ordered = lfc[match(ordered.c1, rownames(lfc)),]
lfc.c2.ordered = lfc[match(ordered.c2, rownames(lfc)),]

n1 = length(c1)
n2 = length(c2)
lfc.c1.ordered.melt = melt(lfc.c1.ordered)
lfc.c1.ordered.melt$Days = c(rep(1,n1), rep(6,n1), rep(12,n1), rep(17,n1), rep(22,n1), rep(27,n1))
colnames(lfc.c1.ordered.melt)[3] = "logFC(asthmatic/healthy)"
lfc.c1.ordered.melt$miRNA = c(rep(seq(1,n1,1),6))

lfc.c2.ordered.melt = melt(lfc.c2.ordered)
lfc.c2.ordered.melt$Days = c(rep(1,n2), rep(6,n2), rep(12,n2), rep(17,n2), rep(22,n2), rep(27,n2))
colnames(lfc.c2.ordered.melt)[3] = "logFC(asthmatic/healthy)"
lfc.c2.ordered.melt$miRNA = c(rep(seq(1,n2,1),6))

col = c("#08FF00",  "#00FF28", 
                 "#00FF58",  "#00FF87",  "#00FFB7", 
                 "#00FFE7",  "#00E7FF", 
                 "#00B7FF",  "#0087FF", "#0058FF", 
                 "#0028FF", "#0800FF")
                 
col.df = data.frame(colors = col,
                    membership = seq(0.55, 0.78, length.out = 12))

col1 = character()
col2 = character()
for (i in 1:length(membership.c1$Cluster1)){
  temp = membership.c1$Cluster1[i]
  col1 = c(col1, col.df$colors[which(abs(col.df$membership-temp)==
                                       min(abs(col.df$membership-temp)))])
}
for (i in 1:length(membership.c2$Cluster2)){
  temp = membership.c2$Cluster2[i]
  col2 = c(col2, col.df$colors[which(abs(col.df$membership-temp)==
                                       min(abs(col.df$membership-temp)))])
}


g.c1 =
  ggplot(lfc.c1.ordered.melt, aes(x = Days, y = `logFC(asthmatic/healthy)`, colour = Var1)) + 
  geom_line() + 
  geom_line(data = fpca.c1 , aes(x=x, y=y), 
            size = 2, color='red', linetype = 2) + 
  scale_colour_manual(name="",values=col1,
                      guide = "none")+
  ggtitle("Upregulated Cluster")+
  # xlab("Day after RV challenge") +
  xlab("") +
  scale_x_continuous(limits = c(0, 30),
                     breaks = seq(0, 30, by = 6))+
  annotate("text", x=0, y=-1.9, label="RV", color="red",
           fontface="bold") +
  annotate("segment", x=0, y=-2.1, xend=0, yend=-2.5,
           col="red", arrow=arrow(length=unit(0.3, "cm")),
           size=1) +
  coord_cartesian(ylim = c(-2.5, 2), clip="off") +
  theme_classic()+
  theme(axis.title = element_text(size = 14, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.5))



g.c2 = ggplot(lfc.c2.ordered.melt, aes(x = Days, y = `logFC(asthmatic/healthy)`, colour = Var1)) + 
  geom_line() + 
  geom_line(data = fpca.c2, aes(x=x, y=y), 
            size = 2, color='red', linetype = 2) + 
  scale_colour_manual(name="",values=col2,
                      guide = "none")+
  ggtitle("Downregulated Cluster")+
  xlab("Day after RV challenge") + 
  scale_x_continuous(limits = c(0, 30),
                     breaks = seq(0, 30, by = 6))+
  annotate("text", x=0, y=-1.9, label="RV", color="red",
           fontface="bold") +
  annotate("segment", x=0, y=-2.1, xend=0, yend=-2.5,
           col="red", arrow=arrow(length=unit(0.3, "cm")),
           size=1) +
  coord_cartesian(ylim = c(-2.5, 2), clip="off") +
  theme_classic()+
  theme(axis.title = element_text(size = 14, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.5))


g.c1.1 = arrangeGrob(g.c1, top = textGrob(expression(bold("A")), x = unit(0.03, "npc"),
                                          y = unit(.7, "npc"), just=c("left","top"),
                                          gp=gpar(col="black", fontsize=20, fontfamily="sans")))

g.c2.1 = arrangeGrob(g.c2, top = textGrob(expression(bold("B")), x = unit(0.03, "npc"),
                                          y = unit(.7, "npc"), just=c("left","top"),
                                          gp=gpar(col="black", fontsize=20, fontfamily="sans")))

# png(filename="../Data/colorkey.png", width = 4, height = 8, units = "in", res = 300)
par(mar = c(2, 2, 15, 15))
image(1,seq(0, 256, length.out = length(col)), 
      t(seq_along(col)), col=col, axes=FALSE,
      xlab = "", ylab = "")
for (i in 1:12){
  text(1.25,seq(0, 256, length.out = 12)[i],
       as.character(format(round(seq(0.78, 0.55, length.out = 12), 2), nsmall = 2)[i]),
       adj = -1, xpd=NA, cex=1.5)}
text(1.5, 296, "Membership", cex = 2, xpd = NA, font = 2)

# dev.off()


colorkey = readPNG("../Data/colorkey.png")
grob = rasterGrob(colorkey)
lay = rbind(c(1,1,1,NA),
            c(1,1,1,NA),
            c(2,2,2,3),
            c(2,2,2,3))

# png("./mfuzz/mfuzz1.png", width = 10, height = 8, units = "in", res = 300)
grid.arrange(grobs = list(g.c1.1, g.c2.1, grob),
             layout_matrix = lay)
# dev.off()





