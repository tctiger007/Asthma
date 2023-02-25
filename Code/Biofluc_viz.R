path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# source("./Biofluc_gss_dataLoader.R")


################# ################# ################# ################# 
######################## check the pvals ###############################
################# ################# ################# ################# 
pvals.deseq2 = read.csv("../Data/deseq2_pval1-289.csv", check.names = FALSE)
colnames(pvals.deseq2)[1] = "miRNAs"
pvals.deseq2.rmwarn = pvals.deseq2[which(pvals.deseq2$warn == "No warn"), ]

n = dim(pvals.deseq2.rmwarn)[1]
colnames(pvals.deseq2.rmwarn)
#0.05/n = 0.0001992032
post_vs_pre_in_asthma = pvals.deseq2.rmwarn$miRNAs[which(pvals.deseq2.rmwarn$post_vs_pre_in_asthma < 0.05/n)]          # 64
post_vs_pre_in_healthy = pvals.deseq2.rmwarn$miRNAs[which(pvals.deseq2.rmwarn$post_vs_pre_in_healthy < 0.05/n)]        # 49

asthma_vs_healthy_at_post = pvals.deseq2.rmwarn$miRNAs[which(pvals.deseq2.rmwarn$asthma_vs_healthy_at_post < 0.05/n)]  # 26
asthma_vs_healthy_at_pre = pvals.deseq2.rmwarn$miRNAs[which(pvals.deseq2.rmwarn$asthma_vs_healthy_at_pre < 0.05/n)]    # 0

pvals.deseq2.rmwarn$miRNAs[which.min(pvals.deseq2.rmwarn$post_vs_pre_in_asthma)]     # "hsa-let-7a-5p"
pvals.deseq2.rmwarn$miRNAs[which.min(pvals.deseq2.rmwarn$post_vs_pre_in_healthy)]    # "hsa-let-7a-5p"
pvals.deseq2.rmwarn$miRNAs[which.min(pvals.deseq2.rmwarn$asthma_vs_healthy_at_post)] # "hsa-miR-122-5p"
pvals.deseq2.rmwarn$miRNAs[which.min(pvals.deseq2.rmwarn$asthma_vs_healthy_at_pre)]  # "hsa-miR-411-3p"

intersect(post_vs_pre_in_healthy, post_vs_pre_in_asthma)  #45
intersect(asthma_vs_healthy_at_pre, asthma_vs_healthy_at_post) # 0 



set.seed(3)
# ### Visulization (give examples to show the following significant miRNAs based on their normalized expression)
# # 1. show post vs pre is significant in healthy group 
miRNA1 = "hsa-miR-122-5p"; group = "H"; nn = 20

count = get(miRNA1, CountsWithGroupInfo.normalized)
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
df = new
df = cbind.data.frame(df, fitted) #, fitted.u95, fitted.l95)
df.group = df[df$Subject_Group == group,]
colnames(df.group)[6] = "fitted.H"
df.group.m = melt(df.group, id.vars = c("Time", "Subject_Group", "Run", "Subject_ID", "InfectionStatus"))

temp = CountsWithGroupInfo.normalized %>% filter(Subject_Group == group) %>%
  select(c("Time", "Subject_Group", "Run", "Subject_ID", "InfectionStatus", all_of(miRNA1)))
temp.m = melt(temp, id.vars = c("Time", "Subject_Group", "Run", "Subject_ID", "InfectionStatus"))

df.group.m = rbind.data.frame(df.group.m, temp.m)
df.group.m$variable = factor(df.group.m$variable,
                             levels = c("hsa-miR-122-5p", "fitted.H"))

spline.df.temp = spline(df.group.m %>% filter(variable == "fitted.H") %>% select(c("Time", "value")))
spline.df = data.frame(spline.df.temp)
spline.df = cbind.data.frame(spline.df, rep("fitted.H", dim(spline.df)[1]))
colnames(spline.df) = c("Time", "value", "variable")


p.group1 =
  ggplot(df.group.m, aes(x = Time, y = value, color = variable)) +
  geom_point() +
  geom_line(data = spline.df, size = 1) +
  xlab("Days") +
  ylab(paste(miRNA1, " \n (Normalized Count)")) +
  # ggtitle("Within-Healthy \n pre- vs. post-RV challenge in healthy subjects") +
  # ggtitle("Group 1: Differential miRNA expression between the virus-challenged \n phase and the control phase in healthy subjects") +
  scale_color_manual("", values = c("hsa-miR-122-5p" = "grey57",
                                    "fitted.H" = "#00ABFD"),
                     guide = guide_legend(override.aes = list(
                       linetype = c("blank", "solid"),
                       shape = c(rep(16, 2))))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-60, 30),
                     breaks = c(seq(-60, 0, by = 10),
                                0, seq(0, 30, by = 5))) +
  scale_y_continuous(limits = c(min(df.group.m$value)*.95, max(df.group.m$value)*1.05),
                     breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = .5),
        axis.title.x = element_text(size = 14, face = "bold",
                                    margin = margin(2,0,0,0, unit = "mm")),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5,
                                  margin = margin(-3,3,2,0, unit = "mm")),
        plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size = 14),
        legend.position = c(1, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(-10,0,0,0),
        legend.box.margin = margin(0,0,0,0),
        legend.box.background = element_rect(colour = "black"))
  
# png(paste0(miRNA1, "Group1.png"),
# height = 8, width = 20, units = "in", res = 300)
p.group1
# dev.off()


# # 2. show post vs pre is significant in asthmatic group 
miRNA2 = "hsa-let-7f-5p"; group = "A"; nn = 20

count = get(miRNA2, CountsWithGroupInfo.normalized)
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
# fitted.u95 = exp(est1$fit + 1.96 * est1$se.fit)
# fitted.l95 = exp(est1$fit - 1.96 * est1$se.fit)
df = new
df = cbind.data.frame(df, fitted) #, fitted.u95, fitted.l95)
df.group = df[df$Subject_Group == group,]
colnames(df.group)[6] = "fitted.A"
df.group.m = melt(df.group, id.vars = c("Time", "Subject_Group", "Run", "Subject_ID", "InfectionStatus"))

temp = CountsWithGroupInfo.normalized %>% filter(Subject_Group == group) %>%
  select(c("Time", "Subject_Group", "Run", "Subject_ID", "InfectionStatus", all_of(miRNA2)))
temp.m = melt(temp, id.vars = c("Time", "Subject_Group", "Run", "Subject_ID", "InfectionStatus"))

df.group.m = rbind.data.frame(df.group.m, temp.m)
df.group.m$variable = factor(df.group.m$variable,
                             levels = c("hsa-let-7f-5p", "fitted.A"))

spline.df.temp = spline(df.group.m %>% filter(variable == "fitted.A") %>% select(c("Time", "value")))
spline.df = data.frame(spline.df.temp)
spline.df = cbind.data.frame(spline.df, rep("fitted.A", dim(spline.df)[1]))
colnames(spline.df) = c("Time", "value", "variable")

p.group2 = ggplot(df.group.m, aes(x = Time, y = value, color = variable)) +
  geom_point() +
  geom_line(data = spline.df, size = 1) +
  xlab("Days") +
  ylab(paste(miRNA2, " \n (Normalized Count)")) +
  # ggtitle("Group 2: Differential miRNA expression between the virus-challenged \n phase and the control phase in asthmatic subjects") +
  # ggtitle("Within-Asthmatic \n pre- vs. post-RV challenge in asthmatic subjects") +
  scale_color_manual("", values = c("hsa-let-7f-5p" = "black",
                                    "fitted.A" = "#FF689F"),
                     guide = guide_legend(override.aes = list(
                       linetype = c("blank", "solid"),
                       shape = c(rep(16, 2))))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-60, 30),
                     breaks = c(seq(-60, 0, by = 10),
                                0, seq(0, 30, by = 5))) +
  scale_y_continuous(limits = c(min(df.group.m$value)*.95, max(df.group.m$value)*1.05),
                     breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = .5),
        axis.title.x = element_text(size = 14, face = "bold",
                                    margin = margin(2,0,0,0, unit = "mm")),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5,
                                  margin = margin(-3,3,2,0, unit = "mm")),
        plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size = 14),
        legend.position = c(1, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(-10,0,0,0),
        legend.box.margin = margin(0,0,0,0),
        legend.box.background = element_rect(colour = "black"))


#
# png(paste0(miRNA2, "Group2.png"),
#     height = 8, width = 20, units = "in", res = 300)
p.group2
# dev.off()




# 3. show asthmatics and healthy controls are not significantly different at pre-RV 
miRNA3 = "hsa-miR-101-3p"; nn = 20

## only plot post-RV period 
count = get(miRNA3, CountsWithGroupInfo.normalized)
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
df = new
df = cbind.data.frame(df, fitted) #, fitted.u95, fitted.l95)
df.m = melt(df, id.vars = c("Time", "Subject_Group", "Run", "Subject_ID", "InfectionStatus"))

temp = CountsWithGroupInfo.normalized %>%
  select(c("Time", "Subject_Group", "Run", "Subject_ID", "InfectionStatus", all_of(miRNA3)))
temp.m = melt(temp, id.vars = c("Time", "Subject_Group", "Run", "Subject_ID", "InfectionStatus"))

df.m = rbind.data.frame(df.m, temp.m)
df.m$variable = paste0(df.m$variable, ".", df.m$Subject_Group)
df.m$variable = factor(df.m$variable,
                       levels = c(paste0(miRNA3, ".H"), 
                                  paste0(miRNA3, ".A"),
                                  "fitted.H", "fitted.A"))

spline.df.tempA = spline(df.m %>% filter(variable == "fitted.A") %>% 
                           select(c("Time", "value")))
spline.dfA = data.frame(spline.df.tempA)
spline.dfA = cbind.data.frame(spline.dfA, rep("fitted.A", dim(spline.dfA)[1]))
colnames(spline.dfA) = c("Time", "value", "variable")

spline.df.tempH = spline(df.m %>% filter(variable == "fitted.H") %>% 
                           select(c("Time", "value")))

spline.dfH = data.frame(spline.df.tempH)
spline.dfH = cbind.data.frame(spline.dfH, rep("fitted.H", dim(spline.dfH)[1]))
colnames(spline.dfH) = c("Time", "value", "variable")

spline.df = rbind.data.frame(spline.dfA, spline.dfH)
# fitHvalue = df.m %>% dplyr::filter(Time == 0 & variable == "fitted.H") %>% dplyr::select(value)
# fitAvalue = df.m %>% dplyr::filter(Time == 0 & variable == "fitted.A") %>% dplyr::select(value)
# 
# addRows = data.frame(Time=c(0,0), 
#                      value=c(fitHvalue,fitAvalue), 
#                      variable=c("fitted.H","fitted.A"))
# spline.df = rbind.data.frame(spline.df, addRows)
p.group3 =  ggplot(df.m %>% filter(Time<=0), aes(x = Time, y = value, color = variable)) + 
  geom_point() +
  geom_line(data = spline.df %>% filter(Time<=0), size = 1) +
  xlab("Days") +
  ylab(paste(miRNA3, " \n (Normalized Count)")) +
  # ggtitle("Between-Unchallenged: \n asthamtic vs. healthy subjects before RV challenge") + 
  # ggtitle("Group 3: Differential miRNA expression between asthmatic \n and healthy subjects at the virus-challenged phase") + 
  scale_color_manual("", values = c("hsa-miR-101-3p.H" = "grey57",
                                    "hsa-miR-101-3p.A" = "black",
                                    "fitted.H" = "#00ABFD",
                                    "fitted.A" = "#FF689F"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(rep("blank",2), rep("solid",2)),
                       shape = c(rep(16, 4))))) +
  # geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(-60, 30),
                       breaks = c(seq(-60, 0, by = 10),
                                  0, seq(0, 30, by = 5))) +
    scale_y_continuous(limits = c(0,2000),
                       breaks = c(seq(0, 2000, by = 500))) +
  theme_classic() +
  coord_cartesian(xlim =c(-60, 0), clip = "on") +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = .5),
        axis.title.x = element_text(size = 14, face = "bold", 
                                    margin = margin(2,0,0,0, unit = "mm")),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, 
                                  margin = margin(-3,3,2,0, unit = "mm")),
        plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size = 14),
        legend.position = c(1, 1),
        legend.justification = c("right", "top"),
        # legend.box.just = "right",
        legend.margin = margin(-15,0,0,0),
        legend.box.margin = margin(0,0,0,0),
        legend.box.background = element_rect(colour = "black"))


# png(paste0(miRNA3, "Group3.png"),
#     height = 8, width = 20, units = "in", res = 300)
p.group3
# dev.off()


# 4. show asthmatics and healthy controls are significantly different at post-RV 
miRNA4 = "hsa-let-7a-5p"; nn = 20

## only plot post-RV period 
count = get(miRNA4, CountsWithGroupInfo.normalized)
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
df = new
df = cbind.data.frame(df, fitted) #, fitted.u95, fitted.l95)
df.m = melt(df, id.vars = c("Time", "Subject_Group", "Run", "Subject_ID", "InfectionStatus"))

temp = CountsWithGroupInfo.normalized %>%
  select(c("Time", "Subject_Group", "Run", "Subject_ID", "InfectionStatus", all_of(miRNA4)))
temp.m = melt(temp, id.vars = c("Time", "Subject_Group", "Run", "Subject_ID", "InfectionStatus"))

df.m = rbind.data.frame(df.m, temp.m)
df.m$variable = paste0(df.m$variable, ".", df.m$Subject_Group)
df.m$variable = factor(df.m$variable,
                       levels = c("hsa-let-7a-5p.H", "hsa-let-7a-5p.A",
                                  "fitted.H", "fitted.A"))

spline.df.tempA = spline(df.m %>% filter(variable == "fitted.A") %>% 
                           select(c("Time", "value")))
spline.dfA = data.frame(spline.df.tempA)
spline.dfA = cbind.data.frame(spline.dfA, rep("fitted.A", dim(spline.dfA)[1]))
colnames(spline.dfA) = c("Time", "value", "variable")

spline.df.tempH = spline(df.m %>% filter(variable == "fitted.H") %>% 
                           select(c("Time", "value")))

spline.dfH = data.frame(spline.df.tempH)
spline.dfH = cbind.data.frame(spline.dfH, rep("fitted.H", dim(spline.dfH)[1]))
colnames(spline.dfH) = c("Time", "value", "variable")

spline.df = rbind.data.frame(spline.dfA, spline.dfH)

p.group4 = ggplot(df.m %>% filter(Time>0), aes(x = Time, y = value, color = variable)) + 
  geom_point() +
  geom_line(data = spline.df %>% filter(Time>0), size = 1) +
  xlab("Days") +
  ylab(paste(miRNA4, " \n (Normalized Count)")) +
  # ggtitle("Between-Challenged: \n asthamtic vs. healthy subjects after RV challenge") + 
  # ggtitle("Group 4: Differential miRNA expression between asthmatic \n and healthy subjects at the virus-challenged phase") + 
  scale_color_manual("", values = c("hsa-let-7a-5p.H" = "grey57",
                                    "hsa-let-7a-5p.A" = "black",
                                    "fitted.H" = "#00ABFD",
                                    "fitted.A" = "#FF689F"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(rep("blank",2), rep("solid",2)),
                       shape = c(rep(16, 4))))) +
  scale_x_continuous(limits = c(-60, 30),
                     breaks = c(seq(-60, 0, by = 10),
                                0, seq(0, 30, by = 5))) +
  scale_y_continuous(limits = c(0,8000),
                     breaks = c(seq(0, 8000, by = 2000))) +
  theme_classic() +
  coord_cartesian(xlim =c(0, 30), clip = "on") +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = .5),
        axis.title.x = element_text(size = 14, face = "bold", 
                                    margin = margin(2,0,0,0, unit = "mm")),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, 
                                  margin = margin(-3,3,2,0, unit = "mm")),
        plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = c(1, 1),
        legend.justification = c("right", "top"),
        legend.margin = margin(-15,0,0,0),
        legend.box.margin = margin(0,0,0,0),
        legend.box.background = element_rect(colour = "black"))


# png(paste0(miRNA4, "Group4.png"),
#     height = 8, width = 20, units = "in", res = 300)
p.group4
# dev.off()


## just for graphical abstract
## don't need real data 
spline.df1 = spline.df
spline.df1$variable = ifelse(spline.df1$variable == "fitted.H", "Healthy", "Asthmatic")
spline.df1$value[118] = 4280
spline.df1$value[117] = 4340


ggplot(data = spline.df1 %>% filter(Time>0), 
       aes(x = Time, y = value, color = variable)) +
  geom_line(size = 2) +
  xlab("Days after rhinovirus challenge") +
  ylab("miRNA4 expression") +
  ggtitle("miRNA") +
  scale_color_manual("", values = c("Healthy" = "#00ABFD",
                                    "Asthmatic" = "#FF689F"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(rep("solid",2)),
                       shape = c(rep(16, 2))))) +
  scale_x_continuous(limits = c(-60, 30),
                     breaks = c(seq(-60, 0, by = 10),
                                0, seq(0, 30, by = 5))) +
  # scale_y_continuous(limits = c(0,8000),
  #                    breaks = c(seq(0, 8000, by = 2000))) +
  theme_classic() +
  coord_cartesian(xlim =c(0, 30), clip = "on") +
  theme(axis.title.y = element_text(size = 30, face = "bold", hjust = .5),
        axis.title.x = element_text(size = 30, face = "bold", 
                                    margin = margin(2,0,0,0, unit = "mm")),
        plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        title = element_text(face = "bold", size = 30),
        plot.title = element_text(hjust = 0.5))

## Venn Diagram
group1_miRNAs = read_excel("../Data/sig_gss.xlsx", sheet = 1)
group2_miRNAs = read_excel("../Data/sig_gss.xlsx", sheet = 2)
group3_miRNAs = read_excel("../Data/sig_gss.xlsx", sheet = 3)
group1_miRNAs = group1_miRNAs$miRNAs
group2_miRNAs = group2_miRNAs$miRNAs
group3_miRNAs = group3_miRNAs$miRNAs
# # 
# # # miRNAs group1 group2 and group3
venn.miRNA = venn.diagram(
  x = list(group1_miRNAs, group2_miRNAs, group3_miRNAs),
  category.names = c("Within-Healthy" , 
                     "Within-Asthmatic", 
                     "Between-Challenged"),
  filename = NULL,
  # './LRT_Pvals_CPMandDESeq2/miRNA_venn.png',
  # output=TRUE,

  # Output features
  imagetype="png" ,
  compression = "lzw",

  # Circles
  lwd = 2,
  col = c("purple", "olivedrab3", "darkorange"),
  fill = c(alpha("purple",0.3), alpha("olivedrab3", 0.3), alpha("darkorange", 0.3)),

  # Numbers
  cex = 1.3,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-10, 10, 0),
  cat.dist = c(0.055, 0.055, 0.035),
  cat.fontfamily = "sans",
  cat.col = c("purple", "olivedrab3", "darkorange"))
#

# 
# png("./miRNA_venn.png", unit = "in", height = 5, width = 10, res = 300)
# grid.arrange(gTree(children = venn.miRNA), top = "Differentially expressed miRNAs")
# dev.off()
# 
#### plot three examples on the same figure
lay = rbind(c(1,1,1,1,2,2,2,2),
            c(1,1,1,1,2,2,2,2),
            c(1,1,1,1,2,2,2,2),
            c(1,1,1,1,2,2,2,2),
            c(1,1,1,1,2,2,2,2),
            c(3,3,3,4,4,4,5,5),
            c(3,3,3,4,4,4,5,5),
            c(3,3,3,4,4,4,5,5),
            c(3,3,3,4,4,4,5,5),
            c(3,3,3,4,4,4,NA,NA))

p.group1.1 = arrangeGrob(p.group1, top = textGrob(expression(bold("A")), x = unit(0.03, "npc"),
                                                  y = unit(.7, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=20, fontfamily="sans")))
p.group2.1 = arrangeGrob(p.group2, top = textGrob(expression(bold("B")), x = unit(0.03, "npc"),
                                                  y = unit(.7, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=20, fontfamily="sans")))

p.group3.1 = arrangeGrob(p.group3, top = textGrob(expression(bold("C")), x = unit(0.03, "npc"),
                                                  y = unit(.7, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=20, fontfamily="sans")))

p.group4.1 = arrangeGrob(p.group4, top = textGrob(expression(bold("D")), x = unit(0.03, "npc"),
                                                  y = unit(.7, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=20, fontfamily="sans")))

venndiagram.1 = arrangeGrob(gTree(children=venn.miRNA), top = textGrob(expression(bold("E")), x = unit(0.03, "npc"),
                                                  y = unit(.7, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=20, fontfamily="sans")))

# png(paste0("./", "sigExamples7.png"),
    height = 15, width = 20, units = "in", res = 300)
grid.arrange(
             p.group3.1, p.group4.1, p.group1.1, p.group2.1,venndiagram.1,
             layout_matrix = lay)
# dev.off()


lay1 = rbind(c(1,1,1,1,2,2,2,2),
             c(1,1,1,1,2,2,2,2),
             c(1,1,1,1,2,2,2,2))
p.group3.2 = arrangeGrob(p.group3, top = textGrob(expression(bold("A")), x = unit(0.03, "npc"),
                                                  y = unit(.7, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=20, fontfamily="sans")))

p.group4.2 = arrangeGrob(p.group4, top = textGrob(expression(bold("B")), x = unit(0.03, "npc"),
                                                  y = unit(.7, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=20, fontfamily="sans")))

# png("figure3.png",
#     height = 8, width = 20, units = "in", res = 300)
grid.arrange(p.group3.2, p.group4.2,             
             layout_matrix = lay1)
# dev.off()


lay2 = rbind(c(1,1,1,1,2,2,2,2,NA,NA),
             c(1,1,1,1,2,2,2,2,3,3),
             c(1,1,1,1,2,2,2,2,NA,NA))
venndiagram.2 = arrangeGrob(gTree(children=venn.miRNA), top = textGrob(expression(bold("C")), x = unit(0.03, "npc"),
                                                                       y = unit(.7, "npc"), just=c("left","top"),
                                                                       gp=gpar(col="black", fontsize=20, fontfamily="sans")))

# png("supple_figure3.png",
#     height = 8, width = 20, units = "in", res = 300)
grid.arrange(p.group1.1, p.group2.1, venndiagram.2,           
             layout_matrix = lay2)
# dev.off()
