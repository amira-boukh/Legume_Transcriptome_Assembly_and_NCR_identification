library(ggplot2)
library(tidyverse)

IndigoRoots_table =  read.table("IndigoRoots_salmon_NCRs.txt",header = F)
IndigoNods_table =  read.table("IndigoNods_salmon_NCRs.txt",header = F)

colnames(IndigoNods_table) = c("Gene_ID","tpm","plant")
colnames(IndigoRoots_table) = c("Gene_ID","tpm","plant") 

IndigoNods_table$logtpm = log10(IndigoNods_table$tpm)
IndigoNods_table$plant = factor(IndigoNods_table$plant)
IndigoNods_table$type = rep("putative NCR",length(IndigoNods_table[,1]))

IndigoRoots_table$logtpm = log10(IndigoRoots_table$tpm)
IndigoRoots_table$plant = factor(IndigoRoots_table$plant)
IndigoRoots_table$type = rep("putative NCR",length(IndigoRoots_table[,1]))

all_df = rbind(IndigoNods_table,IndigoRoots_table)

IndigoNods_combine_df = IndigoNods_table %>%
  filter(plant == "IndigoNods")

IndigoRoots_combine_df = IndigoRoots_table %>%
  filter(plant == "IndigoRoots")

pairwise.wilcox.test(all_df$logtpm,all_df$plant)

ggplot(all_df,aes(x = plant, y= logtpm,fill = plant)) +
  geom_boxplot(position = position_dodge(width = 0.2), width = 0.5)+
  scale_fill_manual(values = c("#FABF7B","#6CB0D6"), name = "Plant", label = c(expression(italic("IndigoNods")), expression(italic("IndigoRoots"))))+
  theme(panel.background = element_rect(
    fill = "white",
    colour = "black",
    linetype = "solid",
  ), legend.position="bottom",axis.text.x = element_blank(),text = element_text(size=20))+
  xlab("NCRs")+
  ylim(-2.5,5)+
  ylab("log10 TPM (Transcripts Per Million)")+
  scale_x_discrete(expand=c(0,1))+
  annotate("text", x = 1.5, y = 5, label = "***", size = 5)+
  annotate("segment", x = 1, xend = 1, y = 4.3, yend = 5)+
  annotate("segment", x = 2, xend = 2, y = 2.5, yend = 5)+
  annotate("segment", x = 1, xend = 2, y = 5, yend = 5)

ggsave("IndigoNCRs_TPM.png",width = 10, height = 7, dpi = 600)



