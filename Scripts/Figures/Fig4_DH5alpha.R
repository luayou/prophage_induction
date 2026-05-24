#Checking with DH5alpha, take an 1day old colony into 5ml of LB + chlora 16h ON, dil 1:50 into 25ml LB +chlor, 
##gre ~2h until 0.4-0.5 OD, collect the cells (5min 12,000g) wash them into new fresh media (with or w/o aTc), monitor OD

rm(list=ls())
library(readxl)
library(tidyverse)
library(ggpubr)
library(ggtext)


#===== Loading files =====

colors_fig <- c("#ff62d9","#628bff","#ffd662")
raw_data <- read_xlsx("PATH_TO_RAW_DATA_PLATE_READER")

mapping_file <- read_xlsx("PATH_TO_MAPPING_FILE")


#==== Parsing the data ====

raw_data$min_time <-  seq(0,(nrow(raw_data)-1)*5,by=5)
info_longer <- raw_data %>% pivot_longer(cols = mapping_file$Well, names_to = "Well", values_to="OD_600") %>%
  select(min_time,Well,OD_600) %>% mutate(h_time=min_time/60)


info_longer <- merge(info_longer, mapping_file)

info_longer <- info_longer[which(info_longer$Well!="E5"),]
info_longer <- info_longer[which(info_longer$Well!="G6"),]

info_longer$Insert <- factor(info_longer$Insert,levels = c("HicA","HicAB","DhiT","DhiTA","Empty", "C(-)"))

info_longer <- info_longer[which(info_longer$aTc=="T"),]

#=== PLOTS ===

H2 <- ggplot(info_longer[is.element(info_longer$Insert,c("HicA","HicAB","Empty"))&info_longer$Isolate_ID==1,],
            aes(h_time,OD_600,colour=Insert)) +# geom_jitter(size=0.01,alpha=0.8) +
  geom_smooth(method = "loess",linewidth = 0.4,alpha=1,fill="gray60") + # stat_smooth() +
  scale_color_manual(values = colors_fig) +
   facet_grid(~aTc, labeller = as_labeller(c("T"="Wilby 2"))) +theme_light() + 
  labs(x="Time [h]",y="*E. coli* DH5&alpha; + pGL001 [OD<sub>600</sub>]")+
  theme(legend.position = 'bottom',axis.title.x = element_markdown(size = 7),axis.text.x = element_text(size = 5,angle = 0),
        axis.title.y = element_markdown(size = 7),axis.text.y = element_text(size=5), 
        strip.background = element_rect(fill="white"), 
        strip.text=element_text(color="black",size=7),legend.title = element_blank(),legend.text = element_text(size=7),
        legend.key.height = unit(3,"mm"))

H2


d2 <- ggplot(info_longer[is.element(info_longer$Insert,c("DhiT","DhiTA","Empty"))&info_longer$Isolate_ID==1,],
             aes(h_time,OD_600,colour=Insert)) + #geom_point(size=0.01,alpha=0.8) + 
  facet_grid(~aTc, labeller = as_labeller(c("T"="Wilby 1")))  +
  geom_smooth(method = "loess", linewidth = 0.4,alpha=1,fill="gray60")+  #stat_smooth() + 
  theme_light() + 
  scale_color_manual(values = colors_fig) + 
  labs(x="Time [h]",y="*E. coli* DH5&alpha; + pGL001 [OD<sub>600</sub>]") +
  theme(legend.position = 'bottom',axis.title.x = element_markdown(size = 7),axis.text.x = element_text(size = 5,angle = 0),
        axis.title.y = element_markdown(size = 7),axis.text.y = element_text(size=5), 
        strip.background = element_rect(fill="white"), 
        strip.text=element_text(color="black",size=7),legend.title = element_blank(),legend.text = element_text(size=7),
        legend.key.height = unit(3,"mm"))

d2
lower_plot <- ggarrange(d2,H2,common.legend = F, nrow = 1, ncol = 2)

ggsave(lower_plot, width = 16.5, height = 7 , units = 'cm', dpi = 320,
       filename = "PATH_TO_FIG4_LOWER_PANEL")
   