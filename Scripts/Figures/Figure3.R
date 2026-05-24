#== Figure 3a ====


#========= Setting up the space =======

rm(list = ls())
invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)),   # Unload add-on packages
                 detach,
                 character.only = TRUE, unload = TRUE))

library(ggplot2)
library(readxl)
library(tidyr)
library(tidyverse)
library(data.table)
library(ggpubr)


#====== Loading files =======

cov_all <- fread("PATH_TO_BP_DEPTH_PHAGES_FROM__ALL_LIBRARIES")
cov_all <- cov_all[!is.na(cov_all$Treatment),]

stats_0 <- read.delim("PATH_TO_MEAN_READS_PER_REGION_PER_LIBRARY")
stats_0 <- stats_0 %>% mutate(RPK=mean_cov/(len_region/1000),perCov=len_cov/len_region*100)
mapping <- read_delim("PATH_TO_MAPPING_FILE")
library_stats <- read.delim2("PATH_TO_TOTAL_NUMBER_OF_READS_PER_LIBRARY",
                             stringsAsFactors = F)

phages_regions <- read.csv("PATH_TO_BED_PHAGE_REGIONS")
phages_regions <- phages_regions[-which(phages_regions$PhageID==""),]

#=== Parsing mapping file ===

mapping$Sample_ID <- gsub("[.]","_",mapping$Sample_ID)

#=== Parsing phage regions ===

phages_regions <- phages_regions[-grep("Cryptic",phages_regions$Phage),]

phages_regions$len <- phages_regions$end- phages_regions$start

#==== Parsing library stats ===

library_stats <- library_stats[grepl("*_1P",library_stats$Sample),]
library_stats$library <- sapply(library_stats$Sample,function(x){strsplit(x,split = "_1P")[[1]][1]})
library_stats$fastqc.total_sequences <- as.numeric(as.character(library_stats$fastqc.total_sequences))
library_stats$Sample_ID <- sapply(library_stats$Sample,function(x){strsplit(x,split = "_S")[[1]][1]})

mapping <- merge(mapping,library_stats,all.x = T)

mapping <- mapping %>% select(Treatment,Replicate,Time_point_h,fastqc.total_sequences) %>%
  mutate(seqs=paste(round(fastqc.total_sequences,1),"M",sep = ""))

ggplot(mapping[is.element(mapping$Treatment,c("YCFA+MMC","YCFA")),],aes(Time_point_h,fastqc.total_sequences,colour=Treatment)) +
  geom_point(size=2)



#=== only using MMC and control ===

cov_all <- cov_all[is.element(cov_all$Treatment,c("YCFA+MMC","YCFA")),]
cov_all <- cov_all[!is.element(cov_all$region,c("contig00002","contig00004")),]

cov_all$Mbp <- cov_all$bp/1000000 
cov_all$region <- factor(cov_all$region,levels = c("Wilby_1" ,"LoVE","Wilby_2","Sombra","Hanky","Shia"))



colors_fig <- c("#72A86B","#A16BA8")

# Prepare a named vector for the new facet names
time_labels <- c("6.5" = "5.5", "10" = "9", "25" = "24") # The time in the mapping file was time since the subculture, adjust to add MMC
time_labels_h <- c("6.5" = "5.5 h", "10" = "9 h", "25" = "24 h") 


#=== Binning ===

binds_width <- 500
flanked_region_bp <- 10000

region=c()
bin_end=c()

for(i in 1:nrow(phages_regions)){
  
  bins_phage <- c()
  bins_by_width <- seq(from=phages_regions$start[i]-flanked_region_bp,
                       to=phages_regions$end[i]+flanked_region_bp,by=binds_width)
 
   if(phages_regions$end[i]%% binds_width ==0){
    
    bins_phage <- bins_by_width[-1]
    
  }else{
    
    bins_phage <- c(bins_by_width[-1],phages_regions$end[i]+flanked_region_bp)
    
  }
  
  region <- c(region,rep(phages_regions$PhageID[i],length(bins_phage)))
  bin_end <- c(bin_end,bins_phage)
}

bin_all_phages <- data.frame(region=region,bin_end=bin_end)
rm(region,bin_end)

samples <- mapping[is.element(mapping$Treatment,c("YCFA+MMC","YCFA")),]

cov_all_bins <- merge(bin_all_phages,samples,all = T)

cov_all_bins$cov_norm_mean_bin <- NA
cov_all_bins$pos_mean_bin <- NA


for(i in 1:nrow(cov_all_bins)){
  
  print(i)
  previous_end <- NA
  
  phage_i <- cov_all_bins$region[i]
  start_phage <- phages_regions$start[which(phages_regions$PhageID==phage_i)] -
    flanked_region_bp
  
  if(cov_all_bins$bin_end[i]==start_phage+binds_width){
    previous_end <- start_phage
  } else{ previous_end <- cov_all_bins$bin_end[i-1] }
  subset_cov_all <- cov_all %>% filter(region == phage_i,
                                       Treatment == cov_all_bins$Treatment[i],
                                       Replicate == cov_all_bins$Replicate[i],
                                       Time_point_h ==cov_all_bins$Time_point_h[i],
                                       bp >= previous_end & bp < cov_all_bins$bin_end[i])
 
  cov_all_bins$cov_norm_mean_bin[i] <- mean(subset_cov_all$cov_norm)
  cov_all_bins$pos_mean_bin[i] <- mean(subset_cov_all$bp)
}

rm(cov_all)
cov_all_bins$Mbp <- cov_all_bins$pos_mean_bin/1000000 
cov_all_bins$region <- factor(cov_all_bins$region,levels = c("Wilby_1" ,"LoVE","Wilby_2","Sombra","Hanky","Shia"))

write.table(cov_all_bins,"PATH_TO_OUT_DEPTH_PER_BINS_TABLE",
            sep = "\t",row.names = F, quote = F)

# === Figure COV per bin ===

p <- ggplot(cov_all_bins,aes(Mbp,cov_norm_mean_bin,colour = Treatment)) + 
  geom_point(size = 0.001,alpha=0.8) + 
  facet_grid(Time_point_h~region,scales = "free_x",space = "free_x",
             labeller = labeller(Time_point_h=time_labels_h)) + theme_bw() + 
  scale_color_manual(values=colors_fig,labels=c("Control","MMC 0.3 ug/mL")) +
  scale_y_continuous(name= "Mapped reads per million of total reads") +
  scale_x_continuous(name= "Bt_806 genomic location [M bp]") +
  theme(legend.position = "bottom",legend.box = "vertical",
        legend.text = element_text(size=7), legend.title = element_text(size=7),
        axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 1,size = 5),
        axis.title.x = element_text(vjust = 0,size = 7),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.background = element_rect(fill="white"), 
        strip.text=element_text(color="black",size=7))


p

ggsave(filename = "PATH_TO_OUTPUT_COVERAGE_PLOT_PER_BIN", plot = p,device = "png")

#qpcrs

library(devtools)
library(scales) # to access break formatting functions
source_gist("524eade46135f6348140") #https://gist.github.com/kdauria/524eade46135f6348140

#install.packages("ggh4x")
library(ggh4x)

#===== Data: CI-like phages  ======

mp_all_love <- read.csv("PATH_TO_MAPPING_FILE_LOVE_QPCRS")
cp_all_love <- read.delim("PATH_TO_CP_VALUES_LOVE_QPCRS",skip = 1,sep = "\t",stringsAsFactors = F)

mp_all_sombra <- read.csv("PATH_TO_MAPPING_FILE_SOMBRA_QPCRS")
cp_all_sombra <- read.delim("PATH_TO_CP_VALUES_SOMBRA_QPCRS",skip = 1,sep = "\t",stringsAsFactors = F)

dna_per_sample <- read_xlsx("PATH_TO_QUBIT_RESULTS")

#==== Total DNA ===

dna_per_sample <- dna_per_sample[which(dna_per_sample$Treatment=="YCFA+MMC" | dna_per_sample$Treatment=="YCFA"),]
dna_per_sample$Jodee <- as.numeric(as.character(dna_per_sample$Jodee))
ggplot(dna_per_sample,
       aes(as.factor(Time_point_h),as.numeric(Jodee),fill= Treatment, shape = as.factor(Replicate))) + 
  geom_bar(stat = "identity",position = "dodge") + facet_grid(~Replicate)

mean_samples <- dna_per_sample %>% group_by(Time_point_h,Treatment) %>% 
  summarise(mean_biological=mean(Jodee),sd=sd(Jodee)) %>% ungroup()


#=== panel 3C ====

h0 <- ggplot(mean_samples,aes(as.factor(Time_point_h),mean_biological,colour = Treatment)) + 
  geom_point(position=position_dodge(width = 0.9),size=1) +
  geom_errorbar(data = mean_samples, aes(as.factor(Time_point_h) , 
                                         y = mean_biological, ymin = mean_biological-sd, ymax = mean_biological+sd), 
                position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.0),width=0.2,alpha=1,linewidth = 0.3) +
  geom_point(data=dna_per_sample,aes(as.factor(Time_point_h),Jodee),alpha=0.4,
             position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.0)) +
  theme_bw() +
  scale_color_manual(values = colors_fig) +
  scale_y_continuous(name=sprintf('Total extracted DNA ng/\u03BCl]')) +
  scale_x_discrete(name = "Time after adding MMC [h]",labels=time_labels) +
  theme(axis.text.x = element_text(size = 6, angle = 0),
        legend.position = "none",axis.title.x = element_text(vjust = 0,size = 7),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7))

h0


dna_log_increase <- dna_per_sample %>% ungroup() %>% select(Replicate,Treatment,Time_point_h,Jodee) %>% 
  group_by(Replicate,Treatment,Time_point_h) %>% pivot_wider(names_from = Treatment,values_from =Jodee) %>% mutate(increase= `YCFA+MMC`/YCFA)

dna_log_increase <- mean_samples %>% ungroup() %>% select(Treatment,Time_point_h,mean_biological) %>% 
  group_by(Treatment,Time_point_h) %>% pivot_wider(names_from = Treatment,values_from =mean_biological) %>% mutate(increase= `YCFA+MMC`/YCFA)


#==== Merging info ===

love <- merge(mp_all_love,cp_all_love,by = "Pos")
sombra <- merge(mp_all_sombra,cp_all_sombra,by = "Pos")

#=== Standards ===

love_standards <- love[which(love$Sample_type=="Standard"),]
love_model <- lm(Cp~Concentration_standard_ng_ul,data = love_standards)

ggplot(love_standards[which(love_standards$Concentration_standard_ng_ul!=1e-8),],aes(log10(Concentration_standard_ng_ul),Cp)) + ggtitle("LoVE") +
  stat_smooth_func(geom = "text",method = "lm",hjust=0,parse=T,xpos = -4,ypos = -5) +
  geom_smooth(method = "lm") + theme_bw() +
  geom_point() #+ facet_wrap(~Primer) 


sombra_standards <- sombra[which(sombra$Sample_type=="Standard"),]
sombra_model <- lm(Cp~Concentration_standard_ng_ul,data = sombra_standards)

ggplot(sombra_standards,aes(log10(Concentration_standard_ng_ul),Cp)) + ggtitle("Sombra") +
  stat_smooth_func(geom = "text",method = "lm",hjust=0,parse=T,xpos = -4,ypos = -5) +
  geom_smooth(method = "lm") + theme_bw() +
  geom_point() #+ facet_wrap(~Primer) 

#=== Using the standard curve ===

love_samples <- love[which(love$Sample_type=="DNA"),]
sombra_samples <- sombra[which(sombra$Sample_type=="DNA"),]

ggplot(love_samples,
       aes(as.factor(Extraction_time_h_after_adding._MMC),Cp)) + geom_point() + facet_grid(~Treatment) + scale_y_log10()

ggplot(sombra_samples,
       aes(as.factor(Extraction_time_h_after_adding._MMC),Cp,shape = as.factor(Technicall_replicate))) + geom_point() + facet_grid(~Treatment) + scale_y_log10()

love_samples$DNA_concentration <- 10^((love_samples$Cp-6.65)/-4.33) 
sombra_samples$DNA_concentration <- 10^((sombra_samples$Cp-8.11)/-3.93) 

ggplot(love_samples,
       aes(as.factor(Extraction_time_h_after_adding._MMC),DNA_concentration,colour = Treatment,shape = as.factor(Technicall_replicate))) + geom_point() + scale_y_log10()

ggplot(sombra_samples,
       aes(as.factor(Extraction_time_h_after_adding._MMC),DNA_concentration,colour = Treatment,shape = as.factor(Technicall_replicate))) + geom_point() + scale_y_log10()

#===== Summary both phages =====

both <- rbind(love_samples,sombra_samples)

#removing unconsistent technical replicates
both <- both[-which(both$Primers=="LoVE"&both$Replicate==2&both$Extraction_time_h_after_adding._MMC==9&both$Technicall_replicate==1),]
both <- both[-which(both$Primers=="LoVE"&both$Replicate==2&both$Extraction_time_h_after_adding._MMC==24&both$Technicall_replicate==2),]

#merging average technical replicates

both_mean_con <- both %>% group_by(Extraction_time_h_after_adding._MMC,Primers,Replicate,Treatment) %>% summarise(mean_technical=mean(DNA_concentration)*10) # the 10 factor reflects the dilution


ggplot(both_mean_con,
       aes(as.factor(Extraction_time_h_after_adding._MMC),mean_technical,colour = Primers)) + 
  geom_point() + scale_y_log10() + facet_grid(Replicate~Treatment)

ggplot(both_mean_con,
       aes(Primers,mean_technical,fill = Treatment)) + 
  geom_bar(stat = "identity",position = "dodge")  + facet_grid(Replicate~as.factor(Extraction_time_h_after_adding._MMC))

ggplot(both_mean_con,
       aes(Treatment,mean_technical,fill = Primers)) + 
  geom_bar(stat = "identity",position = "dodge")  + facet_grid(Replicate~as.factor(Extraction_time_h_after_adding._MMC))

mean_rep <- both_mean_con %>% group_by(Extraction_time_h_after_adding._MMC,Primers,Treatment) %>% summarise(mean_biological=mean(mean_technical),sd=sd(mean_technical)) %>% ungroup()
both_mean_con <- both_mean_con %>% ungroup()


h3 <- ggplot(mean_rep,aes(as.factor(Extraction_time_h_after_adding._MMC),mean_biological,colour = Treatment)) + 
  geom_point(position=position_dodge(width = 0.9),size=1) +
  geom_errorbar(data = mean_rep, aes(as.factor(Extraction_time_h_after_adding._MMC) , y = mean_biological, ymin = mean_biological-sd, ymax = mean_biological+sd),
                position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.0),width=0.2,alpha=1,linewidth = 0.3) +
  geom_point(data=both_mean_con,aes(as.factor(Extraction_time_h_after_adding._MMC),mean_technical),alpha=0.6,size=1,
             position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.0)) +
  facet_grid(~ Primers ,scales = "free_x")+ theme_bw() +
  scale_color_manual(values = colors_fig) +
  scale_y_continuous(breaks = c(1e-3,1e-2, 1e-1),
                     labels = trans_format("log10", math_format(10^.x)),name=sprintf('Phage DNA concentration [ng/\u03BCl]')) +  scale_x_discrete(name="Time after adding MMC [h]") +
  theme(legend.position = "bottom",legend.text = element_text(size=7), legend.title = element_text(size=7),
        axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 1,size = 5),
        axis.title.x = element_text(vjust = 0,size = 7),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.background = element_rect(fill="white"), 
        strip.text=element_text(color="black",size=7))


h3

low_panel <- ggarrange(p,h3,ncol = 2,nrow = 1,widths = c(2,1.1))
low_panel

ggsave(low_panel, width = 17.5, height = 11 , units = 'cm', dpi = 320,
       filename = "PATH_TO_LOW_PANEL_FIG3")

ggsave(h0, width = 6, height = 7 , units = 'cm', dpi = 320,
       filename = "PATH_TO_QUBIT_FIG_3")

