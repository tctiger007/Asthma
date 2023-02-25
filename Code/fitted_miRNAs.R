library(rstudioapi) # only useful in Rstudio
library(gss)
library(readxl)

set.seed(1)
path = dirname(getActiveDocumentContext()$path)
setwd(path)
# 1. show post vs pre is sigsnificant in healthy group: Group 1
# 2. show post vs pre is significant in asthmatic group: Group 2
# 3. show asthmatics and healthy controls are significant different at post-RV: Group 3 
# 4. show asthmatics and healthy controls are significant different at the baseline: Group 4 (no miRNAs)
DE.group1 = read_excel("../Data/sig_gss.xlsx", sheet = 1)
DE.group2 = read_excel("../Data/sig_gss.xlsx", sheet = 2)
DE.group3 = read_excel("../Data/sig_gss.xlsx", sheet = 3)
DE.group1 = DE.group1$miRNAs  # 49
DE.group2 = DE.group2$miRNAs  # 64
DE.group3 = DE.group3$miRNAs  # 26
### DE.group4 0

CountsWithGroupInfo.normalized = read.csv("../Data/CountsWithGroupInfo_normalized_forgss.csv",
                                          check.names = FALSE,
                                          row.names = 1)
CountsWithGroupInfo.normalized = CountsWithGroupInfo.normalized[, -1]
factor.idx = which(colnames(CountsWithGroupInfo.normalized) %in% c("Subject_ID","Subject_Group","Visit","DaysPoint", 
                                                                   "Run", "InfectionStatus", "SampleBins", "Asthma"))
CountsWithGroupInfo.normalized = cbind.data.frame(CountsWithGroupInfo.normalized[,-factor.idx], 
                                                  lapply(CountsWithGroupInfo.normalized[, factor.idx],
                                                         as.factor))
CountsWithGroupInfo.normalized = CountsWithGroupInfo.normalized[,c(1:3, 
                                                                   (dim(CountsWithGroupInfo.normalized)[2]-7):dim(CountsWithGroupInfo.normalized)[2],
                                                                   4: (dim(CountsWithGroupInfo.normalized)[2]-8))]
miRNAs = colnames(CountsWithGroupInfo.normalized)[12:dim(CountsWithGroupInfo.normalized)[2]]

sigs = read.csv("../Data/deseq2_pval1-289.csv", check.names = FALSE,
                row.names = 1)
sigs = sigs[which(sigs$warn == "No warn"), ]
n = dim(sigs)[1]
sigs = sigs[which(sigs$post_vs_pre_in_asthma < 0.05/n | sigs$post_vs_pre_in_healthy < 0.05/n |
                    sigs$asthma_vs_healthy_at_post < 0.05/n | sigs$asthma_vs_healthy_at_pre < 0.05/n),]

DE.sigs = rownames(sigs)
fit.df = data.frame()
nn = 20
for (i in 1:length(DE.sigs)){
  set.seed(3)
  count = get(DE.sigs[i], CountsWithGroupInfo.normalized)
  fit1 = gssanova(count ~ Time + Subject_Group + InfectionStatus + InfectionStatus:Subject_Group + Run,
                  data = CountsWithGroupInfo.normalized, seed = 3,
                  family = "poisson", random = ~1|Subject_ID)
  asthmatics = levels(CountsWithGroupInfo.normalized$Subject_ID)[grepl("A", levels(CountsWithGroupInfo.normalized$Subject_ID))]
  healthy = levels(CountsWithGroupInfo.normalized$Subject_ID)[grepl("H", levels(CountsWithGroupInfo.normalized$Subject_ID))]
  new = data.frame(Time = rep(seq(min(CountsWithGroupInfo.normalized$Time),
                                  max(CountsWithGroupInfo.normalized$Time),
                                  length.out = nn), 2),
                   Subject_Group = as.factor(c(rep("A", nn), rep("H", nn))),
                   Run = as.factor(sample(levels(CountsWithGroupInfo.normalized$Run), 2*nn, replace = TRUE)),
                   Subject_ID = as.factor(c(sample(asthmatics, nn, replace = TRUE),
                                            sample(healthy, nn, replace = TRUE))))
  new$InfectionStatus = factor(ifelse(new$Time<=0, "0", "1"))
  est1 = predict(fit1, new, se = TRUE)
  fitted = exp(est1$fit)
  fit.df = rbind(fit.df, fitted)
}
fit.df = as.data.frame(t(as.matrix(fit.df)))
colnames(fit.df) = DE.sigs
Time = rep(seq(min(CountsWithGroupInfo.normalized$Time),
               max(CountsWithGroupInfo.normalized$Time),
               length.out = nn), 2)
fit.df$Time = Time
rownames(fit.df) = NULL
Subject_Group = as.factor(c(rep("A", nn), rep("H", nn)))
fit.df$Subject_Group = Subject_Group
fit.df$InfectionStatus = factor(ifelse(new$Time<=0, "0", "1"))
fit.df = fit.df[, c((dim(fit.df)[2]-2):dim(fit.df)[2],1:(dim(fit.df)[2]-3))]
write.csv(fit.df, "../Data/fitted_DEmiRNA.csv")