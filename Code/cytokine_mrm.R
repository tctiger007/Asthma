library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(ggpubr)
library(grid)
library(gridExtra)
library(latex2exp)

set.seed(1)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

########################################################
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plot.ft = function(ft = c("IFNg", "IP10", "IL13", "IL17A", "IL6",
                "IL8", "TNFa", "IL33", "IL1b", "IL10",
                "FVC", "FEV1", "PEF",
                "FeNO", "eosinophils", "neutrophils")){
  label = if(ft == "IFNg"){expression(bold(paste("log(IFN-",gamma,") pg/mL")))
  }else if(ft == "IL1b"){expression(bold(paste("log(IL-1",beta,") pg/mL")))
  }else if(ft == "TNFa"){expression(bold(paste("log(TNF-",alpha,") pg/mL")))
  }else if(ft == "FVC"){"log(FVC%)"
  }else if(ft == "FEV1"){expression(bold(paste("log(", FEV[1]*"%", ")")))
  }else if(ft == "PEF"){("log(PEF%)")
  }else if(ft == "eosinophils"){"log(Eosinophil%)"
  }else if(ft == "neutrophils"){"log(Neutrophil%)"
  }else if(ft == "IP10"){"log(IP-10) pg/mL"
  }else if(ft == "IL13"){"log(IL-13) pg/mL"
  }else if(ft == "IL17A"){"log(IL-17A) pg/mL"
  }else if(ft == "IL6"){"log(IL-6) pg/mL"
  }else if(ft == "IL8"){"log(IL-8) pg/mL"
  }else if(ft == "IL33"){"log(IL-33) pg/mL"
  }else if(ft == "IL10"){"log(IL-10) pg/mL"
  }else{"log(FeNO)"}
  
  label1 = if(ft == "IFNg"){TeX("IFN-$\\gamma\\,\\dagger\\ddagger$", bold = TRUE)#expression(bold(paste("IFN-",gamma,""^"\\aleph")))
  }else if(ft == "IL1b"){TeX("IL-1$\\beta\\,\\ddagger\\aleph$", bold = TRUE)   #expression(bold(paste("IL-1",beta,"")))
  }else if(ft == "TNFa"){TeX("TNF-$\\alpha\\,\\ddagger$", bold = TRUE)   #expression(bold(paste("TNF-",alpha,"")))
  }else if(ft == "FVC"){TeX("FVC$\\,\\ddagger$", bold = TRUE)
  }else if(ft == "FEV1"){TeX("FEV$_1\\,\\ddagger$", bold = TRUE) #expression(bold(FEV[1]\\aleph""))
  }else if(ft == "PEF"){TeX("PEF", bold = TRUE)
  }else if(ft == "eosinophils"){TeX("Eosinophil$\\,\\aleph$", bold = TRUE)
  }else if(ft == "neutrophils"){TeX("Neutrolphil$\\,\\dagger\\ddagger\\aleph$", bold = TRUE)
  }else if(ft == "IP10"){TeX("IP-10$\\ddagger$", bold = TRUE)
  }else if(ft == "IL13"){TeX("IL-13$\\,\\ddagger$", bold = TRUE)
  }else if(ft == "IL17A"){TeX("IL-17A$\\,\\dagger\\ddagger\\aleph$", bold = TRUE)#expression(bold(paste("IL-17A"^"\\aleph")))
  }else if(ft == "IL6"){TeX("IL-6$\\,\\dagger\\ddagger$", bold = TRUE)#expression(bold(paste("IL-6"^"\\aleph")))
  }else if(ft == "IL8"){TeX("IL-8$\\,\\ddagger$", bold = TRUE)
  }else if(ft == "IL33"){TeX("IL-33", bold = TRUE)
  }else if(ft == "IL10"){TeX("IL-10$\\,\\ddagger\\aleph$", bold = TRUE)
  }else{TeX("FeNO $\\,\\aleph$", bold = TRUE)}

    
  ft.df = readxl::read_excel(paste0("./", ft, ".xlsx"))
  ft.df$Days = ifelse(ft.df$Time == 1, 0, 
                     ifelse(ft.df$Time == 2, 1,
                            ifelse(ft.df$Time == 3, 3,
                                   ifelse(ft.df$Time == 4, 6, 
                                          ifelse(ft.df$Time == 5, 8, 
                                                 ifelse(ft.df$Time == 6, 10, 13))))))
  ft.df$Group = ifelse(ft.df$Group == "0", "Healthy", "Asthmatic")
  ft.df = na.omit(ft.df)
  ft.fit = ft.df[!duplicated(ft.df$Pred),]
  response = get(paste0("log", ft), ft.df)
  p = ggplot(ft.df, aes(x = Days, 
                    y = get(paste0("log", ft), ft.df), 
                    color = Group)) +
    geom_line(data = ft.fit, 
              aes(x = Days, y = Pred, color = Group)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Group),
                alpha=0.2,
                linetype="dashed",
                color="grey") +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_fill_manual(values = c("#FF689F", "#00ABFD")) +
    scale_colour_manual(values = c("#FF689F", "#00ABFD")) +
    ylab(label) +
    xlab("Day after RV challenge") +
    ggtitle(label1) +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,12,3)) +
    scale_y_continuous(breaks = c(-3,0,5,10), limit = c(-3.5, 10)) +
    theme(axis.title = element_text(face = "bold",size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 13),
          title = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 18),
          legend.position = "bottom")
  }


clinicalFeatures = c(
  "IFNg", "IP10", "IL13", "IL17A", "IL6",
  "IL8", "TNFa", "IL33", "IL1b", "IL10")
g.IFNg = plot.ft(clinicalFeatures[1])
g.IP10= plot.ft(clinicalFeatures[2])
g.IL13 = plot.ft(clinicalFeatures[3])
g.IL17A = plot.ft(clinicalFeatures[4])
g.IL6 = plot.ft(clinicalFeatures[5])
g.IL8 = plot.ft(clinicalFeatures[6])
g.TNFa = plot.ft(clinicalFeatures[7])
g.IL33 = plot.ft(clinicalFeatures[8])
g.IL1b = plot.ft(clinicalFeatures[9])
g.IL10 = plot.ft(clinicalFeatures[10])

legend <- get_legend(g.IFNg)

# grid.newpage()
h <- grobHeight(legend)
w <- grobWidth(legend)
title <- textGrob("Day after RV challenge", y=unit(0.5,"npc") + 0.5*h, 
                  vjust=-.5, gp=gpar(fontsize=18, fontface="bold", fontfamily="sans"))
gt <- gTree(children=gList(legend, title))

# png("../../../../Manuscript/BioFlucFinal/Revision_Viruses_2/Figures/Figure_2_cytokine.png",
#     width = 15, height = 8, units = "in", res = 300)
p.ggarrange = grid.arrange(arrangeGrob(g.IFNg + theme(legend.position = "none", axis.title.x = element_blank()), 
          g.IL6 + theme(legend.position = "none", axis.title.x = element_blank()), 
          g.IL17A + theme(legend.position = "none", axis.title.x = element_blank()), 
          g.IP10 + theme(legend.position = "none", axis.title.x = element_blank()), 
          g.IL1b + theme(legend.position = "none", axis.title.x = element_blank()), 
          g.IL13 + theme(legend.position = "none", axis.title.x = element_blank()),
          g.IL10 + theme(legend.position = "none", axis.title.x = element_blank()), 
          g.IL8 + theme(legend.position = "none", axis.title.x = element_blank()), 
          g.TNFa + theme(legend.position = "none", axis.title.x = element_blank()),
          g.IL33+ theme(legend.position = "none", axis.title.x = element_blank()),
          gt,
          layout_matrix = rbind(c(1,2,3,4,5),
                                c(1,2,3,4,5),
                                c(1,2,3,4,5),
                                c(1,2,3,4,5),
                                c(1,2,3,4,5),
                                c(1,2,3,4,5),
                                c(6,7,8,9,10),
                                c(6,7,8,9,10),
                                c(6,7,8,9,10),
                                c(6,7,8,9,10),
                                c(6,7,8,9,10),
                                c(6,7,8,9,10),
                                c(NA,11,11,11,NA),
                                c(NA,11,11,11,NA))))
# dev.off()




