library(ggplot2)
library(pROC)
library(tidyverse)
library(knitr)
library(ggpubr)
library(MASS)
library(gridExtra)
library(pheatmap)
library(ggVennDiagram)

# Suppress all warning messages
options(warn = -1)
# Prevent creation of Rplots.pdf
options(device = function() pdf(NULL))

# Create Figures/manuscripts directory if it doesn't exist
dir.create("./Figures/manuscripts/", recursive = TRUE, showWarnings = FALSE)

plot.hits<-function(hits, FDR=0.05, FC=1.5, main, 
                    xlim=c(0, 10), ylim=c(-2, 2),
                    colUp="#B31B21", colDown="#1465AC", 
                    xlab="logCPM", ylab="Log(Tumor/Normal)", 
                    FCline=F, ...){
  hits.flt = hits$FDR<FDR
  
  up.flt=hits.flt & hits[, "FC"]>0
  dwn.flt=hits.flt & hits[, "FC"]<0
  plot(hits[, "logCPM"], hits[, "logFC"], 
       pch=20, cex=0.3, col="darkgray", cex.axis=1.8, cex.lab=1.8, font.lab=2,
       main=main, cex.main=1.5,
       ylab=ylab, xlab=xlab, ylim=ylim, xlim=xlim)
  
  points(hits[up.flt, "logCPM"], hits[up.flt, "logFC"],pch=20, cex=0.5, col=colUp)
  points(hits[dwn.flt, "logCPM"], hits[dwn.flt, "logFC"], pch=20, cex=0.5, col=colDown)
  if ( FCline ){
    abline(h=log2(FC), lwd=2,  lty=2)
    abline(h=-log2(FC), lwd=2, lty=2)
  }
  text(9, 1.5, labels=format(sum(up.flt), format="d", big.mark=","), col=colUp, cex=2, font=2)
  text(9, -1.5, labels=format(sum(dwn.flt), format="d", big.mark=","), col=colDown, cex=2, font=2)
}

cohort.colors = c("#B31B21", "#1465AC")

##################################################
# Figure 1A
##################################################
seq_qc_path = "./data/Figure_1_sample_qc.tsv"
df = read.delim(file=seq_qc_path, header=T)

cancer.names= c("Breast", "Colorectal", "Lung", "Ovarian", "Pancreatic", "All") 
names(cancer.names)=c("breast", "colon", "lung", "ovary", "pancreas", "all") 


df=df[with(df, order(primary_diagnosis)),]
df = df %>% dplyr::rename (tissue=primary_diagnosis,
                    Cohort = cohort)

df$Cohort=gsub("Cancer", "Tumor", df$Cohort)
df$Cohort=gsub("Control", "Normal", df$Cohort)

p1a <- ggboxplot(df %>% filter(!is.na(hmc_percent)), x = "tissue", y = "hmc_percent",
               color = "Cohort", palette = rev(cohort.colors),
               add = "jitter", ylab="5hmC mass fraction", xlab="Tissue",
               legend = "right")+
  theme(text = element_text(size = 12, face="bold"))+
  stat_compare_means(aes(group = Cohort), label = "p.signif")

# Save Figure 1A
ggsave("./Figures/manuscripts/Figure1A.pdf", p1a, width=8, height=6)

#############################################################
# Figure 1C
#############################################################
enrichment_path = "./data/Figure_1_gfeature_enrichment.tsv"
df = read.delim(file=enrichment_path, header=T)
df$feature=factor(df$feature,
                  ordered=T,
                  levels=sort(unique(df$feature)))

df = df %>% dplyr::rename (tissue=primary_diagnosis,
                    Cohort = cohort)

df$Cohort=gsub("Cancer", "Tumor", df$Cohort)
df$Cohort=gsub("Control", "Normal", df$Cohort)

names(cohort.colors)=c("Tumor", "Normal")
df$Cohort=factor(df$Cohort,levels=c("Normal", "Tumor"))

plist=list()
for ( g in sort(unique(df$tissue))) {
  
  plist[[g]]<- ggboxplot(subset(df, tissue==g, drop=T),
                         x = "feature", y = "log2Enrich", outlier.shape = NA,
                         color = "Cohort", palette = cohort.colors,
                         xlab="Genomic Feature", ylab="Enrichment (log2)",
                         font.label = list(size = 4, color = "black"),
                         legend = "right") +
    theme(text = element_text(size = 8, face="bold")) +coord_flip()+
    stat_compare_means(aes(group = Cohort), label = "p.signif", hide.ns=T,
                       label.x=3)+ggtitle(g)+theme(
                         plot.title = element_text(hjust = 0.5),
                         legend.position = "none")
}
# create a plot to store legend
tmp<-ggplot(subset(df, tissue==g, drop=T),
            aes(x = feature, y = log2Enrich, col=Cohort)) +geom_boxplot(outlier.shape = NA)+
  scale_color_manual(values=rev(cohort.colors))+
  labs(color = "Cohort")+theme_bw()+
  theme(legend.text=element_text(size=8, face="bold"), 
        legend.title = element_text(face="bold"))
# make legend a separate plot
legend <- ggpubr::get_legend(tmp)
plist[["last"]]=as_ggplot(legend)
grid.arrange(grobs=plist, ncol=3)

########################################################
# Figure 2C
########################################################
drg_per_tissue = "./data/Figure_2_DRG_per_tissue.Rdata"
drg_all_tissue = "./data/Figure_2_DRG_allvsall.Rdata"
load(drg_per_tissue)
load(drg_all_tissue)

pdf("./Figures/manuscripts/Figure2C.pdf", width=10, height=8)
venn.list=lapply(tissue.specific.drg, function(x) rownames(x)[x$FDR<0.05 & x$FC>0])
ggVennDiagram::ggVennDiagram(venn.list, label_alpha = 0,label = "count") +
  scale_fill_gradient(low = "white", high = "white", guide = "none") +
  scale_color_manual(values = rep("black", 5))
dev.off()

########################################################
# Figure 2D
########################################################
pdf("./Figures/manuscripts/Figure2D.pdf", width=10, height=8)
venn.list=lapply(tissue.specific.drg, function(x) rownames(x)[x$FDR<0.05 & x$FC<0])
ggVennDiagram::ggVennDiagram(venn.list, label_alpha = 0,label = "count") +
  scale_fill_gradient(low = "white", high = "white", guide = "none") +
  scale_color_manual(values = rep("black", 5))
dev.off()

#######################################################
# Figure 3B
######################################################
stage_drg_file = "./data/Figure_3_stage_DRG.Rdata"
load(stage_drg_file)

pdf("./Figures/manuscripts/Figure3B.pdf", width=12/1.3, height=8/1.5)

par(mfrow=c(2,3))
par(mar=c(3.8, 4.5, 1.5, 1))
par(oma=c(1.8, 1.8, 1, 1))
par(mgp=c(2.8, 1, 0))  # Reduce spacing between plots while keeping margins

for ( i in c("", ".breast", ".colon", ".lung", ".ovary", ".pancreas")){
  txt=gsub("^\\.", "", ifelse(i!='', i, 'all'))
  xlab = if (i %in% c(".ovary", ".pancreas", ".lung")) {
    expression(log[2]~(CPM))
  } else {
    ""
  }
  ylab = if (i %in% c("", ".lung")) {
    expression(log[2]~(late/early~FC))
  } else {
    ""
  }
plot.hits(stage.early.vs.late.edgR.res[[paste0("III.IVvsI.II", i)]],
            main=txt, xlab=xlab, ylab=ylab)
}
dev.off()
