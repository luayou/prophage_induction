#Sup_fig5

#=== Organizing the environment ===

rm(list=ls())
library(data.table)
library(tidyverse)
library(gggenes)
library(ggpubr)
library(ggrepel)

#==== Loading files =====

a1_genome <- read.delim("PATH_TO_CC00806_ANNOTATIONS")

load("PATH_TO_FCOUNTS_SENSE_R_OBJECT") #This was produced using Rsubread and load fc
sense <- fc 

mapping_file <- read.csv("PATH_TO_MAPPING_FILE",
                         stringsAsFactors = F)
mapping_file <- mapping_file[which(mapping_file$Dowstream=="RNAseq"),]

parsing_bam_id <- data.frame(library=sense$targets,
                             sampleID=sapply(sense$targets,function(x){strsplit(x,split = "_")[[1]][2]}))

mapping_file <- merge(mapping_file,parsing_bam_id,all.x = T)
mapping_file$treatment <- ifelse(mapping_file$treatment=="No mitC","control","mitC")
mapping_file$treatment <- factor(mapping_file$treatment, levels = c("control","mitC"))

mapping_file$prefix_depth <- sapply(mapping_file$library,function(x){strsplit(x,split = "_.bam")[[1]][1]})


#==== Getting depth for all region fwd

counts_sense <-  data.frame(library=sense$counts,strand=sense$annotation)
sense_long <- counts_sense %>% pivot_longer(cols = starts_with("library"),names_to = "library", names_prefix = "library." , values_to = "count" )
sense_long <- merge(sense_long,mapping_file,all.x = T)
sense_long$rpk <- (sense_long$count)/(sense_long$strand.Length/1000)
sense_long <- sense_long %>% group_by(library) %>% mutate (tpm = rpk/(sum(rpk)/1000000))

sense_long <- sense_long %>% ungroup() %>% group_by(strand.GeneID)  %>% mutate(zscore_tpm = (log2(tpm+1) - median(log2(tpm+1)))/sd(log2(tpm+1)) ) 
# z-score based on https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-021-02936-w and https://www.biostars.org/p/107519/ 


ci_like_genes_complete <- sense_long[is.element(sense_long$strand.GeneID,paste("contig00001",c(4904,4905,4906,4907,4908,4909,4910),sep = "_")),]

col_treatments <- c("#72A86B","#A16BA8")
names(col_treatments) <- c("control","mitC")


#ci_like_genes <- sense_long[is.element(sense_long$strand.GeneID,paste("contig00001",c(4689),sep = "_")),]
new_labels <- paste("gp",c(4904,4905,4906,4907,4908,4909,4910),sep = "_")
names(new_labels) <- paste("contig00001",c(4904,4905,4906,4907,4908,4909,4910),sep = "_")


panel_a <- ggplot(ci_like_genes_complete,
                  aes(as.factor(timepointh),asinh(tpm),fill=treatment)) +
  geom_boxplot(linewidth = 0.3) + #geom_hline(yintercept = asinh(23),linetype = "dashed") + 
  facet_wrap(~strand.GeneID,labeller = as_labeller(new_labels))+ scale_fill_manual(values = col_treatments,labels=c("Control","MMC")) +
  scale_y_continuous(name="Asinh(TPM)")+ #,limits = c(3,9),breaks = seq(4,8,2)) +
  scale_x_discrete(name="Time after adding MMC [h]")  +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 1,size = 6),
        legend.position = "bottom",
        axis.ticks.x = element_line(linewidth = 0.3),
        axis.ticks.y = element_line(linewidth = 0.3),
        axis.title.x = element_text(vjust = 0,size = 7),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.background = element_rect(fill="white"), 
        strip.text=element_text(color="black",size=7))


panel_a

to_print <- ggarrange(h,panel_a, ncol = 1, heights = c(1,4))

to_print
ggsave(to_print, width = 18, height = 15, units = 'cm', dpi = 320,
       filename = "PATH_TO_SUP_FIG_5")


