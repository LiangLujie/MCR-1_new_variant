library(ggplot2)
library(ggrepel)
library(dplyr)

# Read in original data
data<-read.delim("M6_vs_WT_MCR-1.txt")
data$log2fc<-log2(data$fc)

data$reg<-ifelse(data$log2fc > 1 & data$pval < 0.05, "Upregulated",
                 ifelse(data$log2fc < -1 & data$pval < 0.05, "Downregulated","nosig"))

# Read in original data
ref<-read.delim("annotated_cluster.txt")
ref<-ref[,c(1,8)]
colnames(ref)[1]="id"

data<-data %>%
  left_join(ref,by="id")

label=c("ldtD","dacC","ampH")

data$label<-ifelse(data$symbol %in% label,data$symbol, NA)

data$result<-ifelse(data$reg == "Upregulated" & !data$symbol %in% label,"up",
                      ifelse(data$reg == "Downregulated" & !data$symbol %in% label,"down",
                             ifelse(data$symbol %in% label,"interest","none")))

# Plot
ggplot(data,aes(x=log2fc,y= -log10(pval)))+
  geom_point(size=4,aes(colour=result))+
  geom_hline(yintercept=-log10(0.05),linetype=4,color="gray",size=1.1) +
  geom_vline(xintercept=c(-1,1),linetype=4,color="gray",size=1.1) +
  scale_color_manual(values=c("blue","chartreuse3","gray","red"))+
  #geom_text_repel(data=data,aes(label=label),color= "chartreuse3", max.overlaps = Inf,size=6)+
  xlab("Log2 Fold Change")+
  ylab("-Log10 (p-value)")+
  theme_classic() + theme(panel.grid = element_blank(),legend.position = "none",
                          axis.text = element_text(size = 20),
                          axis.title = element_text(size = 15),
                          axis.ticks.length = unit(0.2,'cm'),
                          axis.ticks = element_line(size = 2),
                          axis.line = element_line(size = 2))+
  xlim(-4,4)

ggsave("volcano_m6_vs_wt.pdf",device = "pdf",width = 10,height = 8)
