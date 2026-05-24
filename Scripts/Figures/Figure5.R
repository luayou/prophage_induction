#other qpcrs

rm(list=ls())
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(readxl)
library(devtools)
library(scales)
library(ggh4x)

source_gist("524eade46135f6348140") #https://gist.github.com/kdauria/524eade46135f6348140

#===== Data: CI-like phages  ======


summary_table <- read.delim("PATH_TO_PROPHAGE_INFO_LINKED_ISOLATE_CI",
                            stringsAsFactors = F)

pomma <- read.delim("PATH_TO_CP_VALUES_ALL_POMMA",
                    stringsAsFactors = F)

pomma$group_100 <- "B"

pomma$group_100[which(pomma$Isolate=="588")] <- "B"

love <- read.delim("PATH_TO_CP_VALUES_ALL_LOVE",
                    stringsAsFactors = F)


love$group_100 <- "A"

d_761 <- read.delim("PATH_TO_CP_VALUES_CC00761",
                stringsAsFactors = F)


d_761$group_100 <- "C"


ef_881_967 <- read.delim("PATH_TO_CP_VALUES_CC00881ANDCC00967",
                         stringsAsFactors = F)


ef_881_967$group_100 <- "D"
ef_881_967$group_100[which(ef_881_967$Isolate=="967")] <- "E"


all <- rbind(love,pomma,d_761,ef_881_967)


all$isolate <- sprintf("%05d", as.numeric(as.character(all$Isolate)))
all$isolate <- paste("CC",all$isolate,sep = "")

all$isolate <- factor(all$isolate,levels = c("CC00806", "CC00070", "CC01040" ,"CC01534",
                                             "CC00139","CC00340","CC01404","CC01409", "CC00588",
                                             "CC00761","CC00881","CC00967"))

#=== Getting info of other phages and taxa of these isolates ====

isolates_fig_5 <- summary_table[is.element(summary_table$isolate,all$isolate),]

isolates_fig_5$names <- ifelse(isolates_fig_5$cluster_induced=="Yes",isolates_fig_5$induced_name,isolates_fig_5$Cluster_ID_95ani_85cov)
isolates_fig_5$names <- gsub("Cluster","VC",isolates_fig_5$names)
isolates_fig_5$names
isolates_fig_5$names <- gsub("^.*phage* ","",isolates_fig_5$names)
isolates_fig_5$names[which(isolates_fig_5$names == "p00")] <- "Hanky"


group <- all %>% select(isolate,group_100) %>% distinct() 

n_isolates_per_cluster <-  isolates_fig_5 %>% group_by(names) %>% tally() %>% arrange(desc(n))
ci_phages <- c("LoVEphage", "Pomma","VC_217","VC_156","VC_126")
order_phages <- c(ci_phages,n_isolates_per_cluster$names[!is.element(n_isolates_per_cluster$names,ci_phages)])

isolates_fig_5 <- merge(isolates_fig_5,group,all.x = T)

isolates_fig_5$isolate <- factor(isolates_fig_5$isolate,levels = c("CC00806", "CC00070", "CC01040" ,"CC01534",
                                                                       "CC00139","CC00340","CC01404","CC01409", "CC00588",
                                                                       "CC00761","CC00881","CC00967"))


isolates_fig_5$names <- factor(isolates_fig_5$names,levels =order_phages)

color_panel_top <- c("#9e9e9efd","#a8706b")

h <- ggplot(isolates_fig_5,aes(isolate,names,fill = as.factor(n_complete_ci))) +
  geom_tile() + facet_grid(~isolate,scales = "free_x") + theme_bw() +
  scale_fill_manual(values  = color_panel_top) + 
  ylab("Viral cluster")+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill="white"), strip.text=element_text(color="black",size=7),
        panel.grid.major.x = element_blank(),panel.grid.major.y = element_line(linetype = "dashed",linewidth = 0.3),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7))



h

colors_fig <- c("#72A86B","#A16BA8")
#colors_fig <- c("#ffd662","#628bff")

min_v <- min(all$DNA_concentration, na.rm=T)
max_v <- max(all$DNA_concentration, na.rm=T)

p <- ggplot(all,aes(Treatment,DNA_concentration,fill=Treatment)) + geom_boxplot() + theme_bw() +
  scale_fill_manual(values  = colors_fig)  + facet_grid(~isolate) +# + facet_nested(~group_100 + isolate) +
  scale_x_discrete(name="")+
  geom_point(size=0.5) +
  scale_y_log10(limits= c(min_v,max_v),breaks = c(1e-6,1e-5,1e-4,1e-3,1e-2, 1e-1,1),
                                                      labels = trans_format("log10", math_format(10^.x)),name="DNA concentration of virion \n encoding a complete CI-like [ng/ul]") +
  theme(axis.text.x = element_blank(),legend.position = "none",
        axis.ticks.x = element_blank(),
        #strip.background =element_blank(),strip.text = element_blank())
        strip.background = element_rect(fill="white",colour = "white"), strip.text=element_text(color="black",size=7),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7))
p

f <- ggarrange(h,p,ncol = 1,align = "v" )
f

ggsave(filename = "PATH_TO_OUTPUT_FIG5", plot = f ,device = "svg",
       width = 18, height = 14, units = 'cm')

#=== Stats ====

stats <- data.frame(isolate=unique(all$isolate),p_value=NA,shapiro_pvalue=NA)

for(i in 1:nrow(stats)){
  
  control_i <- all[which(all$isolate==stats$isolate[i] & all$Treatment =="Control"),] 
  mmc_i <- all[which(all$isolate==stats$isolate[i] & all$Treatment =="MMC"),]
  if(nrow(control_i) !=  nrow(mmc_i)){
    next
  }
  
  d <- control_i$DNA_concentration - mmc_i$DNA_concentration
  # Shapiro-Wilk normality test for the differences
  d_sh <- shapiro.test(d) # => p-value = 0.91
  stats$shapiro_pvalue[i] <- d_sh$p.value
  
  m <- t.test(control_i$DNA_concentration,mmc_i$DNA_concentration,paired = T, alternative = "less")
  stats$p_value[i] <- m$p.value
}

stats$fdr <- p.adjust(stats$p_value,method = "fdr") 



