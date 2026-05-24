#Merging second and third attempt to show the three induction batches

rm(list=ls())

library(readxl)
library(tidyverse)
library(ggplot2)

library(dplyr)


#===== Loading files =====


raw_gc_25isolates <- read.delim("PATH_TO_RAW_DATA_PLATE_READER_ROUND1",
                                skip = 11,sep = ",")

mapping_file_25 <- read.csv("PATH_TO_MAPPING_FILE1")

raw_gc_47isolates <- read.delim("PATH_TO_RAW_DATA_PLATE_READER_ROUND2",
                                skip = 11,sep = ",")

mapping_file_47 <- read_xlsx("PATH_TO_MAPPING_FILE2")


#==== Parsing the data ====

info_longer_25 <- raw_gc_25isolates %>% pivot_longer(cols = mapping_file_25$Well, names_to = "Well", values_to="OD_600")
info_longer_25 <- info_longer_25 %>% select(Well,Duration..Hours.,Duration..Minutes.,OD_600)

info_longer_25 <- merge(info_longer_25, mapping_file_25)

info_longer_47 <- raw_gc_47isolates %>% pivot_longer(cols = mapping_file_47$Cell, names_to = "Cell", values_to="OD_600")
info_longer_47 <- info_longer_47 %>% select(Cell,Duration..Hours.,Duration..Minutes.,OD_600)

info_longer_47 <- merge(info_longer_47, mapping_file_47)

info_longer_47 <- info_longer_47[!is.element(info_longer_47$Isolate,info_longer_25$Isolate),]

#==== Merging ROUNDS ====

colnames(info_longer_47)[1] <- "Well"
colnames(info_longer_25)[8] <- "MMC_ug_ml"
info_longer_25$MMC_ug_ml[which(info_longer_25$MMC_ug_ml=="wo_MMC")] <- "0"
info_longer <- rbind(info_longer_25,info_longer_47)


colors_fig <- c("#72A86B","#d95f02","#A16BA8")

to_plot_all$Isolate <- factor(to_plot_all$Isolate,levels = c("CC00806", "CC00070", "CC01040" ,"CC01534",
                                                                   "CC00139","CC00340","CC01404","CC01409", "CC00588",
                                                                   "CC00761","CC00881","CC00967"))


p <- ggplot(to_plot_all,aes(Duration..Hours.,OD_600, colour = as.factor(MMC_ug_ml) )) + 
  geom_jitter() + 
  scale_x_continuous(breaks=seq(0,40,4)) + xlim(0,40)  +
  scale_color_manual(values = colors_fig)+#, names= c(0,0.03,0.3),name="MMC [\u03BCg/mL]") +
  facet_wrap(Isolate~.) + theme_bw()  +
  labs(x="Time after adding MMC [h]",colour="MMC [\u03BCg/mL]") + 
  theme(axis.text.x = element_text(size = 6),legend.position = "bottom",
        axis.title.x =  element_text(size = 7),
        strip.background = element_rect(fill="white",colour = "white"), strip.text=element_text(color="black",size=7),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7))


ggsave(p, width = 15, height = 13 , units = 'cm', dpi = 320,
       filename = "PATH_TO_SUP_FIG6")
