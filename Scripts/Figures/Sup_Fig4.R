#Percentage reads per region
rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyverse)

#=== Loading and parsing files ===

regions <- read.csv("PATH_TO_COORDINATES_GENOMIC_REGIONS")
regions <- regions[,1:2]


mapping_file <- read.delim("PATH_TO_MAPPING_FILE",
                               stringsAsFactors = F,skip = 1)
mapping_file$time <- sapply(mapping_file$sample_alias,function(x){strsplit(x,split = "_")[[1]][2]})
mapping_file$treatment <- sapply(mapping_file$sample_alias,function(x){strsplit(x,split = "_")[[1]][3]})
mapping_file$treatment <- factor(mapping_file$treatment,levels = c("No mitC","0.3ug/ml mitC"))
mapping_file$replicate <- sapply(mapping_file$sample_alias,function(x){strsplit(x,split = "_")[[1]][4]})

mapping_file <- mapping_file %>% select(sample_title,time,treatment,replicate)

reads_per_region <- read.delim("PATH_TO_READCOUNT_PER_GENOMIC_REGION_PER_LIBRARY",
                               sep=";",header = F)

colnames(reads_per_region) <- c("library","region","reads")

reads_per_region <- merge(reads_per_region,regions,all.x = T)
reads_per_region$sample_title <- sapply(reads_per_region$library,function(x){strsplit(x,split = "_R1")[[1]][1]})
reads_per_region <- merge(reads_per_region,mapping_file,all.x = T)

reads_per_region <- reads_per_region %>% group_by(sample_title,source,time,treatment,replicate) %>% summarise(total_reads=sum(reads))
reads_per_region <- reads_per_region %>% group_by(sample_title) %>% mutate(per=total_reads/sum(total_reads)*100)
reads_per_region %>% group_by(sample_title) %>% summarise(sum=sum(per))

#== PLOT ===

p <- ggplot(reads_per_region,aes(time, per, colour=source)) + 
  geom_boxplot(whisker.linewidth = 0.1, staple.linewidth = 0.3,median.linewidth = 0.3,box.linewidth = 0.3) + 
  facet_grid(treatment~.) + 
  theme_light() + geom_hline(yintercept = 1,linetype="dashed",colour="gray") +
  theme_bw() + labs(colour="Genomic region") +
  scale_y_sqrt(name=sprintf('Percentage of reads mapping per region')) +
  scale_x_discrete(name = "Time after adding MMC [h]",labels=c(0,1.25,2.25,3.25,4.25)) +
  theme(axis.text.x = element_text(size = 6, angle = 0),
        axis.title.x = element_text(vjust = 0,size = 7),
        legend.title = element_text(size=7),
        legend.text = element_text(size = 6),
        #legend.position = "bottom",
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7))
p

ggsave(p, width = 8, height = 7 , units = 'cm', dpi = 320,
       filename = "PATH_TO_SUPFIG_4A")
