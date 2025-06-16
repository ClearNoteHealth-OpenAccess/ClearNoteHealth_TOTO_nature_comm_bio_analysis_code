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

get.coords.for.ggplot <- function(roc) {
  df <- coords(roc, "all", transpose = FALSE)
  return(df[rev(seq(nrow(df))),])
}

f_logit<-function(p) { log(p/(1-p))}

f_compare_two_SF_model_aucs<-
  function(model1, model2, model1_name, model2_name,
           plot_title=NA, 
           base_font_size=20, 
           paired, 
           alternative="two.sided") {
    
    roc1 = roc(predictor=model1$outerCV_by_sample$Lambda_Min_Test_Response, 
                response=model1$outerCV_by_sample$Label, 
                direction="<", plot=F) 
    roc2 = roc(predictor=model2$outerCV_by_sample$Lambda_Min_Test_Response,
                response=model2$outerCV_by_sample$Label,
                direction="<", plot=F)
    
    roc1$name = model1_name
    roc2$name = model2_name
    
    corr_label<-ifelse(paired==TRUE,"Correlated","")
    
    theme_set(theme_gray(base_size = base_font_size))
    
    if (is.null(roc1$name)) 
      roc1$name<-"Model_1"
    if (is.null(roc2$name))
      roc2$name<-"Model_2"
    if (is.na(plot_title))
      plot_title<-""
    roc_curves_to_plot<-rbind(data.frame("Model"=paste0(roc1$name,"\nauROC=",round(roc1$auc,4)),get.coords.for.ggplot(roc1)),
                              data.frame("Model"=paste0(roc2$name,"\nauROC=",round(roc2$auc,4)),get.coords.for.ggplot(roc2)))
    delong_results<-roc.test(roc1,roc2,
                             method="delong",
                             paired=paired,
                             alternative = alternative)
    # plot both features
    ggplot_obj<- ggplot(roc_curves_to_plot) +
      geom_line(aes(x=1-specificity, y=sensitivity,color=Model), size=1.2) + 
      ylab("Sensitivity") + xlab("1-Specificity") +  
      annotate("text", label = paste0("Two-sided DeLong\np-value: ",signif(delong_results$p.value,2)), x = .75, size=base_font_size/4, y = .25, colour = "black") +
      ggtitle(plot_title)+
      theme_bw()+ theme(axis.line = element_line(color='black'),
                        plot.background = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        axis.text=element_text(size=base_font_size*0.8),
                        axis.title=element_text(size=base_font_size),
                        axis.ticks = element_line(size = 1.2),
                        axis.ticks.length = unit(0.2, "cm"),
                        legend.text=element_text(size=base_font_size*0.8),
                        legend.title=element_text(size=base_font_size*0.8))
    return(ggplot_obj)
  }

plot.one.metagene<-function(dd.ave, i, xlab="Genebody", ylab="5mC Enrichment"){
  d1=subset(dd.ave, tissue==i, drop=T)
  plot(as.numeric(d1[1,3:ncol(dd.ave)]), 
       col=tissue.colors[i], type="l",
       ylab=ylab, ylim=c(-1.2, 0.8),
       xlab=xlab, 
       xaxt='n', yaxt='n', main=i, lwd=2, cex.lab=1.5)
  axis(1, at = c(0, 80), las=1, labels=c("-2K", "TES"), cex=2)
  axis(1, at = c(20, 100), las=1, labels=c("TSS", "+2K"), cex=2)
  axis(2, at = c( -1, -0.5, 0, 0.5, 1, 1.5), las=1, labels=c( -1, -0.5, 0, 0.5, 1, 1.5), cex=2)
  lines(as.numeric(d1[2, 3:ncol(dd.ave) ]), 
        col=tissue.colors[i], 
        type="l", lty=2, lwd=2)
  legend("bottomright",
         legend=c("Tumor", "Normal"),
         col=tissue.colors[i],
         lty=c(1, 2), lwd=2, bty="n") 
}

##################################################
# SFigure 1B
##################################################
seq_qc_path = "./data/Figure_1_sample_qc.tsv"
df = read.delim(file=seq_qc_path, header=T)

cancer.names= c("Breast", "Colorectal", "Lung", "Ovarian", "Pancreatic", "All") 
names(cancer.names)=c("breast", "colon", "lung", "ovary", "pancreas", "all") 

cancer.colors = c('#1f77b4','#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#babab9')
names(cancer.colors) = c("breast", "colon", "lung", "ovary", "pancreas", "all")

df=df[with(df, order(primary_diagnosis)),]
df = df %>% dplyr::rename (tissue=primary_diagnosis,
                    Cohort = cohort)

df$Cohort=gsub("Cancer", "Tumor", df$Cohort)
df$Cohort=gsub("Control", "Normal", df$Cohort)

df.sub =  subset(df, !is.na(age) & !is.na(hmc_percent) & Cohort=='Tumor', drop=T)

# scatter plot 
p_s1b1 <- ggscatter(df.sub, 
          x="age", y="hmc_percent", col="tissue",
          add = "reg.line", xlab="Age", ylab="5hmC Mass Fraction",
          conf.int = T, 
          palette = cancer.colors,
          add.params = list(color = "blue", fill = "lightgray"),
          width = 6, height = 4)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.y = 0.13, size=6) +
  stat_regline_equation(label.y = 0.14, size=6) +
  geom_abline(slope=1, intercept = 0, linetype="dashed", color="red") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
# Save SFigure 1B scatter plot
ggsave("./Figures/manuscripts/SFigure1B.pdf", p_s1b1, width=8, height=6)

##########################################
# SFigure 1C
##########################################
tissue.colors = c('#1f77b4','#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#babab9')
names(tissue.colors)=c("breast", "colon", "lung", "ovary", "pancreas")

metagene_file = "./data/SFigure_1C_metagene.csv"
dd.ave = read.csv(file=metagene_file, check.names = F)[,-1]
colnames(dd.ave)=gsub("^X", "", colnames(dd.ave))

pdf("./Figures/manuscripts/SFigure1C.pdf", width=12, height=3)
par(mfrow=c(1,5))
# Increase left margin to prevent y-label cropping, adjust other margins
par(mar=c(5.1, 4.2, 2.1, 1.1), oma=c(2, 3, 2, 2))
plot.one.metagene(dd.ave, "breast", xlab="")
plot.one.metagene(dd.ave, "colon", "","")
plot.one.metagene(dd.ave, "lung", ylab="")
plot.one.metagene(dd.ave, "ovary","","")
plot.one.metagene(dd.ave, "pancreas","","")
dev.off()

######################################################
# SFigure 1E
#####################################################
c8_gsea = "./data/SFigure_1E_gsea_c8.tsv"
c8_tissue=read.delim(file=c8_gsea,
                     header=T, stringsAsFactors = F)

s1e <- ggplot(c8_tissue, aes(x=subtype, y=mean_fpkm))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.2, size=1, alpha=0.7)+  # Add individual data points
  stat_compare_means(label.y=1.05)+
  xlab("")+ylab("Mean FPKM")+
  theme(axis.text=element_text(size=12, face="bold"),
        strip.text = element_text(size = 12, face="bold"))+
  ylim(0,1.1)+facet_wrap(~pathway, ncol=5)
ggsave("./Figures/manuscripts/SFigure1E.pdf", s1e, width=12, height=3)

###############################################
# SFigure 4A
###############################################

cancer_score_file = "./data/SFigure_4A_cancer_score.tsv"
cancer.info = read.delim(file=cancer_score_file, header=T, stringsAsFactors = F)
s4a <- ggplot(cancer.info, aes(x=stage_current, y=cancer.score))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.1)+
  xlab("Stage")+stat_compare_means(label.y=1.05)

ggsave("./Figures/manuscripts/SFigure4A.pdf", s4a, width=8, height=3)
 
###############################################
# SFigure 4B
###############################################
output_file = "./data/SFigure_4B_models.rds"
figure4sb.models = readRDS(output_file)

non_dhrm_auc = 
  roc(predictor=figure4sb.models[["non_dhmr_model"]]$outerCV_by_sample$Lambda_Min_Test_Response, 
            response=figure4sb.models[["non_dhmr_model"]]$outerCV_by_sample$Label, 
            direction="<", plot=F) 
hyper_dhmr_auc=
   roc(predictor=figure4sb.models[["hyper_dhmr_model"]]$outerCV_by_sample$Lambda_Min_Test_Response, 
          response=figure4sb.models[["hyper_dhmr_model"]]$outerCV_by_sample$Label, 
          direction="<", plot=F)
hypo_dhmr_auc =
  roc(predictor=figure4sb.models[["hypo_dhmr_model"]]$outerCV_by_sample$Lambda_Min_Test_Response,
      response=figure4sb.models[["hypo_dhmr_model"]]$outerCV_by_sample$Label,
      direction="<", plot=F)

p1=f_compare_two_SF_model_aucs(figure4sb.models[["hyper_dhmr_model"]],
                               figure4sb.models[["non_dhmr_model"]], 
                               "Hyper-DhmR",
                               "Non-DhmR",
                               paired=T)
  
p2=f_compare_two_SF_model_aucs(figure4sb.models[["hypo_dhmr_model"]],
                               figure4sb.models[["non_dhmr_model"]],
                               "Hypo-DhmR",
                               "Non-DhmR",
                               paired=T)

combined_plot <- grid.arrange(p1, p2, ncol=2)
ggsave("./Figures/manuscripts/SFigure4B.pdf", combined_plot, width=14, height=4.5)

################################################
# SFigure S4C
################################################
tissue.colors = c('#1f77b4','#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#babab9')

t.color=tissue.colors[c(1,2,2,3,4,5)]
names(t.color)=c("breast", "colon", "colon", "lung", "ovary", "pancreas")
cohort.colors = c("#B31B21", "#1465AC" )
names(cohort.colors)=c("Cancer", "Control")

dd=read.csv(file="./data/SFigure_4C_tm_scores.csv")
rownames(dd)=dd[,1]
dd=dd[,-1]
dd$tissue=factor(dd$tissue, ordered = T,
                 levels=c("breast", "colon","lung", "ovary", "pancreas"))
s4c <- pheatmap(t(dd[, 1:8]),
         cluster_col=F,
         cluster_row=F,
         show_rownames=T,
         show_colnames = F,
         annotation_col=
           data.frame(row.names=rownames(dd),
                      cohort=dd$cohort,
                      tissue=dd$tissue),
         annotation_colors = list(cohort=cohort.colors,
                                  tissue=t.color), fontsize_col=5)
ggsave("./Figures/manuscripts/SFigure4C.pdf", s4c, width=12, height=5)

#####################################
# SFigure 5C
#####################################
dd = read.csv(file="./data/SFigure_5C_tumortissue_toto_tf.csv")
s5c <- ggplot(dd,
       aes(x=correct.pred, y=tumor_fraction))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(col=classification), width=0.1)+
  facet_wrap(~class)
s5c <- s5c + stat_compare_means(
  method = "wilcox.test",
  data = dd %>% filter(class!='ovary.Control'),
  label.y = 0.9,  # Position the text annotation lower in the plot
  label.x = 1  # Position the text annotation lower in the plot
)
ggsave("./Figures/manuscripts/SFigure5C.pdf", s5c, width=12, height=8)
