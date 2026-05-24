# Figure 2a

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

#=== Loading files

annotation_phages <- read_xlsx("PATH_TO_ANNOTATIONS_PHAGES")
annotation_phages$Gene <- sapply(annotation_phages$Gene,function(x){strsplit(x,"AMR_A1_")[[1]][2]})
annotation_phages <- annotation_phages[which(annotation_phages$Gene!="contig00001_2443"),] # this gene doesn't seem to be part of Wilby_2

annotation_phages$Manual_category_wrap <- str_wrap(annotation_phages$Manual_category,20)
unique(annotation_phages$Manual_category_wrap )
annotation_phages$Manual_category_wrap <- factor(annotation_phages$Manual_category_wrap,levels=str_wrap(c("Hypothetical","CI-like protein",
                                                                                                  "moron, auxiliary metabolic gene and host takeover",
                                                                                                  "integration and excision",
                                                                                                  "tRNA" ,
                                                                                                  "5S rRNA",
                                                                                                  "DNA, RNA and nucleotide metabolism",
                                                                                                  "other",
                                                                                                  "connector",
                                                                                                  "tail",
                                                                                                  "head and packaging" ,
                                                                                                  "lysis",
                                                                                                  "transcription regulation"),20))
levels(annotation_phages$Manual_category_wrap)

#=== Colors

df_colors <- data.frame(categories=factor(unique(annotation_phages$Manual_category_wrap),
                                          levels=str_wrap(c("Hypothetical","CI-like protein",
                                          "moron, auxiliary metabolic gene and host takeover",
                                          "integration and excision",
                                          "tRNA" ,
                                          "5S rRNA",
                                          "DNA, RNA and nucleotide metabolism",
                                          "other",
                                          "connector",
                                          "tail",
                                          "head and packaging" ,
                                          "lysis",
                                          "transcription regulation"),20)), colors = NA)

df_colors <- df_colors[order(df_colors$categories),]
                        
df_colors$colors <- c("#bdbdbd","#e31a1c","#ff7f00","#fb9a99",
                      "#fdbf6f","#ffff99","#1f78b4", "#a6cee3",
                      "#b2df8a","#33a02c","#cab2d6","#6a3d9a","#b15928")

#==== Plot

annotation_phages$orientation <- ifelse(annotation_phages$strandedness==1,T,F)

annotation_phages$label <- annotation_phages$Manual
annotation_phages$label[which(annotation_phages$label=="Hypothetical")] <- ""

min_pos_per_phage <- annotation_phages %>% group_by(Phage) %>% summarise(min_loc=min(start_position))

annotation_phages <- merge(annotation_phages,min_pos_per_phage,all.x = T)



annotation_phages$relative_start_position_kb <- (annotation_phages$start_position - annotation_phages$min_loc)/1000
annotation_phages$relative_end_position_kb <- (annotation_phages$end_position - annotation_phages$min_loc)/1000


annotation_phages$PhageID <- factor(annotation_phages$Phage,
                                  levels =c("LoVE","Sombra","Shia","Hanky",
                                            "Wilby_1","Wilby_2"))

phage_labels <- c("LoVE"="LoVE\n~72 K bp\n100%",
                  "Sombra"="Sombra \n~65 K bp\n98.1%",
                  "Shia"= "Shia \n~58 K bp\n98.6%",
                  "Hanky"= "Hanky\n~43 K bp\n100%",
                  "Wilby_1" = "Wilby_1\n~45 k bp\n96.5%",
                  "Wilby_2"= "Wilby_2\n~46 k bp\n99.16%")

length(table(annotation_phages$Manual_category))

p <- ggplot(annotation_phages, aes(xmin = relative_start_position_kb, xmax = relative_end_position_kb, y=PhageID, 
                                   forward = orientation,fill=Manual_category_wrap)) +
  geom_gene_arrow(arrowhead_width = unit(1, "mm"), color="grey30",linewidth = 0.05) + 
  scale_y_discrete(expand = expansion(add = c(16,1)),labels=phage_labels)  +
  scale_fill_manual(values = df_colors$colors)+ 
  facet_wrap(. ~ PhageID, scales = "free_y", ncol = 1) +
  theme_genes()  %+replace% theme(legend.position = "bottom",
                                  legend.text = element_text(size=6),
                                  legend.title = element_blank(),
                                  #legend.key.size = unit(1, "cm"),
                                  axis.text.x = element_text(size = 7), 
                                  axis.text.y = element_text(size = 7),
                                  axis.title.x = element_text(size=7),
                                  axis.ticks.y = element_blank(),
                                  #axis.ticks.x = element_blank(),
                                  plot.margin=unit(x=c(0.5,0,0,0),units="cm"),
                                  panel.grid.major = element_blank(), 
                                  panel.grid.minor = element_blank()) +
  xlab("Prophage genomic location [K bp]") + ylab("") +  coord_cartesian(clip = "off")  +
  geom_text_repel(data=annotation_phages %>% mutate(start2 = (relative_start_position_kb+relative_end_position_kb)/2), 
                  aes(x=start2, label=str_wrap(label,30)),  
                        max.overlaps = Inf,  size = 5 / .pt,
                  angle= 90,nudge_y =-9,
                  segment.linetype = "longdash",
                  segment.size = 0.05) #+ coord_cartesian(clip = "off") 

p

ggsave(p, width = 18, height = 22.5, units = 'cm', dpi = 320,
       filename = "PATH_TO_SUP_FIG2") 
            

