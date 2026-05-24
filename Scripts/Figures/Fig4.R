## Figure 4A

#========== Setting up the workspace ========

rm(list = ls())
invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)),   # Unload add-on packages
                 detach,
                 character.only = TRUE, unload = TRUE))

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gggenes)
library(ggpubr)
library(ggrepel)
library(pheatmap)
library(readxl)

#=== Loading files ===

annotation_phages <- read_xlsx("PATH_TO_PHAGES_ANNOTATION_TABLE")
annotation_phages$Gene <- sapply(annotation_phages$Gene,function(x){strsplit(x,"AMR_A1_")[[1]][2]})

genes_higher10tpm <- read.delim("PATH_TO_PROPHAGE_GENES_HIGHER_10TPM_IN_5_TIMEPOINTS_ACROSS_ALL_SAMPLES_WITHOUT_SIGNIFICANT_DIFF")
length(unique(genes_higher10tpm$strand.GeneID))

annotation_const_allsamples_alltimes <- read.delim("PATH_TO_CURATED_ANNOTATION_OF_CONSTITUTIVE_EXPRESSED_GENES")

#=== Figure genes consist

ggplot(genes_higher10tpm,aes(asinh(tpm),strand.GeneID,colour = treatment)) + geom_boxplot() +geom_vline(xintercept = asinh(10)) + 
  facet_grid(region~timepointh,scales = "free_y",space = "free") + theme_light() #+ theme(axis.text.x = element_blank()) 


#==== Adding_higher10 ==

annotation_phages$higher10in5times <- ifelse(annotation_phages$Gene%in% genes_higher10tpm$strand.GeneID,T,F)
table(annotation_phages$higher10in5times)

#==== Plot ===

annotation_phages$orientation <- ifelse(annotation_phages$strandedness==1,T,F)

annotation_phages <- merge(annotation_phages,annotation_const_allsamples_alltimes,all.x = T)

annotation_phages$PhageID <- factor(annotation_phages$Phage,
                                    levels = c("LoVE","Sombra","Shia","Hanky",
                                               "Wilby_1","Wilby_2"))

annotation_phages <- annotation_phages[which(annotation_phages$Gene!="contig00001_2443"),]

min_pos_per_phage <- annotation_phages %>% group_by(PhageID) %>% summarise(min_loc=min(start_position))

annotation_phages <- merge(annotation_phages,min_pos_per_phage,all.x = T)



annotation_phages$relative_start_position_kb <- (annotation_phages$start_position - annotation_phages$min_loc)/1000
annotation_phages$relative_end_position_kb <- (annotation_phages$end_position - annotation_phages$min_loc)/1000


unique(annotation_phages$Color_category)

annotation_phages$Color_category[is.na(annotation_phages$Color_category)] <- "No_cons"
annotation_phages$Color_category <- factor(annotation_phages$Color_category,levels =unique(annotation_phages$Color_category) )
str(annotation_phages$Color_category)
colors <- c("white","#EC821C","#4d4d4dff","#09BDE0","#1F96FF","#d7191c")#c("#DEDDDE","#dd1c77")#"#B2DF8A")

phage_labels <- c("LoVE"="LoVE\n(~72 K bp)","Sombra"="Sombra \n(~65 K bp)","Shia"= "Shia \n(~58 K bp)",
                  "Hanky"= "Hanky\n(~43 K bp)",
                  "Wilby_1" = "Wilby_1\n(~45 k bp)","Wilby_2"= "Wilby_2\n(~46 k bp)")

max(annotation_phages$relative_end_position_kb)

p <- ggplot(annotation_phages, aes(xmin = relative_start_position_kb, xmax = relative_end_position_kb, y=PhageID, 
                                   forward = orientation,fill=Color_category)) +
  geom_gene_arrow(arrowhead_width = unit(1, "mm"),color="grey30",linewidth = 0.05) + 
  scale_fill_manual(values = colors) +
  scale_y_discrete(expand = expansion(add = c(6,1)),labels=phage_labels) +
  facet_wrap( .~PhageID, scales = "free_y", ncol = 1) +
  xlim(0,71.8) +
  theme_genes()  %+replace% theme(legend.position = "bottom",
                                  #legend.margin = margin(0,0,0,0),
                                  legend.justification = "center",
                                  legend.location = "plot",
                                  legend.text = element_text(size=7),
                                  legend.title = element_text(size=7),
                                  axis.text.x = element_text(size = 5), 
                                  axis.text.y = element_text(size = 7),
                                  axis.title.x = element_text(size=7),
                                  axis.ticks.y = element_blank(),
                                  panel.grid.major = element_blank(), 
                                  panel.grid.minor = element_blank()) +
  xlab("Prophage genomic location [K bp]") + ylab("") +
  geom_text_repel(data=annotation_phages %>% mutate(start2 = (relative_start_position_kb+relative_end_position_kb)/2), 
                  aes(x=start2, label = Region_label), size = 6 / .pt,
                  angle= 0,nudge_y = -2,
                  segment.linetype = "dashed",
                  segment.size = 0.1) + coord_cartesian(clip = "off") 

p

ggsave(p, width = 18, height = 10, units = 'cm', dpi = 320,
       filename = "PATH_TO_FIG4_UPPER_PANEL")
