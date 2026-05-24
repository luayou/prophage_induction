#only bacteroidota predicted phages

rm(list=ls())
library(tidyr)
library(tidyverse)

#=== Loading files ===

summary_table <- read.delim("PATH_TO_PROPHAGE_REGIONS_LINKED_TO_ISOLATE_CI_AND_PREVIOUS_INDUCTIONS_SUP_TABLE",
                          stringsAsFactors = F)

bacteroidota_isolates <- read.delim("PATH_TO_GTDBTK_SUMMARY_SUP_TABLE", 
                                    stringsAsFactors = F)

#=== Getting basic info ====

length(unique(summary_table$isolate)) # 305 isolates with selected phages
round(305/335*100)
length(unique(summary_table$Cluster_ID_95ani_85cov)) # 248 VCs

poly <- summary_table %>% group_by(isolate) %>% tally()
length(poly$isolate[which(poly$n>1)]) #240
round(240/305*100) #~79% of the lysogens are poly
round(240/335*100) #72% of the total bacteroidales collection

induced <- summary_table[(which(!is.na(summary_table$Induced))),]
nrow(induced) # 73 of the selected phages were induced by Sof
clusters_induced <- summary_table[is.element(summary_table$Cluster_ID_95ani_85cov,induced$Cluster_ID_95ani_85cov),]
nrow(clusters_induced) # 310 of the selected phages belong to cluster of at least one member was induced by Dahlman
length(unique(clusters_induced$Cluster_ID_95ani_85cov)) # those 312 phages belong to 22 VCs
length(unique(induced$Cluster_ID_95ani_85cov))


#==== Getting taxa info =====

isolates <- bacteroidota_isolates %>% separate(col=classification,into= c("domain","phylum","class","order","family","genus","species"),
                                               sep = ";")
isolates$phylum <- sapply(isolates$phylum,function(x){strsplit(x,split = "__")[[1]][2]})
isolates$phylum <- sapply(isolates$phylum,function(x){strsplit(x,split = "_")[[1]][1]})

print(isolates$user_genome[which(isolates$phylum!="Bacteroidota")])
isolates <- isolates[which(isolates$phylum=="Bacteroidota"),]

isolates$genus <- sapply(isolates$genus,function(x){strsplit(x,split = "__")[[1]][2]})
isolates$family <- sapply(isolates$family,function(x){strsplit(x,split = "__")[[1]][2]})
isolates$species <- sapply(isolates$species,function(x){strsplit(x,split = "__")[[1]][2]})
table(isolates$genus,useNA = "always")
table(isolates$family,useNA = "always")
table(isolates$order,useNA = "always")

isolates$genus <- sapply(isolates$genus,function(x){strsplit(x,split = "_")[[1]][1]})
isolates$order <- sapply(isolates$order,function(x){strsplit(x,split = "__")[[1]][2]})
isolates$class <- sapply(isolates$class,function(x){strsplit(x,split = "__")[[1]][2]})

isolates <- isolates[,c(1,3,4,5,6,7,8,9)]

#=== 1. Range from genus  ====

unique(isolates$species)
table(isolates$genus)

color_genus <- data.frame(color=NA,genus=unique(isolates$genus))
color_genus$color[which(color_genus$genus=="Alistipes")] <- "#8c2d04"#"#ffffd4"#"#000000" #https://davidmathlogic.com/colorblind/#%23000000-%23E69F00-%2356B4E9-%23009E73-%23F0E442-%230072B2-%23D55E00-%23CC79A7
color_genus$color[which(color_genus$genus=="Bacteroides")] <-"#ffffd4"# "#8c2d04" ##0072B2"
color_genus$color[which(color_genus$genus=="Barnesiella")] <- "#cc4c02"#"#fec44f"##56B4E9"
color_genus$color[which(color_genus$genus=="Coprobacter")] <- "#fee391"#"#E69F00" 
color_genus$color[which(color_genus$genus=="Parabacteroides")] <- "#fe9929"# "#009E73"
color_genus$color[which(color_genus$genus=="Phocaeicola")] <- "#fec44f"#"#cc4c02"##D55E00" # Maybe change this one, for this one "#CC79A7"
color_genus$color[which(color_genus$genus=="Prevotella")] <- "#ec7014"##F0E442"


range_to_itol <- data.frame(isolate=isolates$user_genome,range="range",genus=isolates$genus)
range_to_itol <- merge(range_to_itol,color_genus)
range_to_itol <- range_to_itol %>% relocate(isolate,range,color,genus)

write.table(range_to_itol,"PATH_TO_OUTPUT_RANGE_TABLE",
            sep = ",",row.names = F, col.names = F, quote = F)


#== n phages with complete CI

length(unique(summary_table$Prediction[which(!is.na(summary_table$Prediction))])) #902 prophages
length(unique(summary_table$isolate[which(!is.na(summary_table$Prediction))])) #305 lysogens
305/335*100

length(which(!is.na(summary_table$Prediction))) #902 predicted prophaegs 

length(which(summary_table$n_complete_ci>0)) # 197
length(which(summary_table$n_ntd_ci>0)) # 429

total_phages <- summary_table %>% filter(!is.na(summary_table$Prediction)) %>% 
  group_by(isolate) %>% summarise(n=n(),ci_complete=length(which(n_complete_ci>=1))) 
total_phages$isolate[which(total_phages$isolate=="CC00806")] <- "CC00806T"

isolates_with_complete_ci <- total_phages[which(total_phages$ci_complete>0),]
table(isolates_with_complete_ci$n)

total_phages_0 <- data.frame(isolate=setdiff(isolates$user_genome,total_phages$isolate),n=0,ci_complete=0)

total_phages <- rbind(total_phages,total_phages_0)
length(total_phages$isolate[which(total_phages$n>1)]) # 240 polylysogens
240/305*100

total_phages$no_complete_ci <- total_phages$n - total_phages$ci_complete
total_phages <- total_phages %>% select(isolate,ci_complete,no_complete_ci)


#write.table(total_phages,"PATH_TO_OUTPUT_MULTIBAR_COUNTING",
 #           sep = ",",row.names = F, col.names = F, quote = F)


only_complete_ci <- total_phages %>% select(isolate,ci_complete)
write.table(only_complete_ci,"PATH_TO_OUTPUT_SINGLE_BAR_COUNTING_COMPLETE",
            sep = ",",row.names = F, col.names = F, quote = F)


only_no_complete_ci <- total_phages %>% select(isolate,no_complete_ci)
write.table(only_no_complete_ci,"PATH_TO_OUTPUT_SINGLE_BAR_COUNTING_NON_COMPLETE",
            sep = ",",row.names = F, col.names = F, quote = F)

#=== Figure 1c ====


phages_complete_ci <- summary_table %>% filter(n_complete_ci>0 & !is.na(Prediction))

clusters_n_phages_ci <- phages_complete_ci %>% group_by(Cluster_ID_95ani_85cov,genus,cluster_induced,induced_name) %>% tally()
clusters_n_phages_ci$id <- ifelse(!is.na(clusters_n_phages_ci$induced_name),clusters_n_phages_ci$induced_name,clusters_n_phages_ci$Cluster_ID_95ani_85cov)
clusters_n_phages_ci$id <- sub("Cluster_","VC_",clusters_n_phages_ci$id)
clusters_n_phages_ci$id <- sapply(clusters_n_phages_ci$id,function(x){if(!startsWith(x,"VC")){strsplit(x,split = "phage ")[[1]][2]}else{x}})

p <- ggplot(clusters_n_phages_ci,aes(reorder(id,n),n,fill = genus)) + geom_bar(stat = "identity",position = "stack",colour="black",linewidth = 0.15) + 
  scale_y_continuous(name = "N isolates carrying CI complete prophages") + scale_x_discrete(name = "Viral clusters with phages encoding complete CI") +
  coord_flip() + scale_fill_manual(values = c("#ffffd4","#cc4c02", "#fe9929","#fec44f"),name="Genus") + theme_bw() + theme(legend.position = 'bottom')

r <- ggplot(clusters_n_phages_ci,aes(reorder(id,-n),n,fill = genus)) + geom_bar(stat = "identity",position = "stack",colour="black",linewidth = 0.15) + 
  scale_y_continuous(name = str_wrap("N isolates carrying prophages with complete CI",30)) + scale_x_discrete(name = "Viral clusters with prophages with complete CI") +
  #coord_flip() + 
  scale_fill_manual(values = c("#ffffd4","#cc4c02", "#fe9929","#fec44f"),name="Genus") + theme_light() + 
  theme(legend.position = 'none',axis.title.x = element_text(size=7),axis.text.x = element_text(size = 5,angle = 90),
        axis.title.y = element_text(size=7),axis.text.y = element_text(size=5))#, legend.text = element_text(size=6), legend.title = element_text(size=8))

r
ggsave(p, width = 20, height = 20, units = 'cm',
       filename = "PATH_TO_FIG_1C")
 
ggsave(r, width = 12, height = 5, units = 'cm',
       filename = "PATH_TO_FIG_1C_NOLEGEND")
        

#=== paragraph paper

table(summary_table$n_complete_ci[which(!is.na(summary_table$Prediction))]) # 196 phages with complete CI
length(unique(summary_table$Cluster_ID_95ani_85cov[which(!is.na(summary_table$Prediction) & summary_table$n_complete_ci>0)])) #61 VC

(unique(summary_table$induced_name[which(!is.na(summary_table$Prediction) & summary_table$n_complete_ci>0)])) #9 induced VC

length(unique(summary_table$isolate[which(!is.na(summary_table$Prediction) & summary_table$n_complete_ci>0)])) #154 isolates encode a CI complete phage
round(154/305*100,1)

sum(clusters_n_phages_ci$n) #197 phages with complete ci
sum(clusters_n_phages_ci$n[which(clusters_n_phages_ci$id=="Pomma"|clusters_n_phages_ci$id=="LoVEphage")]) # 60 isolates has either Pomma or Love

60/198*100

love_complete_ci <- summary_table$isolate[which(summary_table$induced_name=="Bacteroides phage LoVEphage" &
                                                  summary_table$n_complete_ci>0)]
love_complete_ci <- c(love_complete_ci,"CC00806T")

love_info <- total_phages[is.element(total_phages$isolate,love_complete_ci),]
table(love_info$ci_complete)

pomma_complete_ci <- summary_table$isolate[which(summary_table$induced_name=="Bacteroides phage Pomma" &
                                                  summary_table$n_complete_ci>0)]
pomma_info <- total_phages[is.element(total_phages$isolate,pomma_complete_ci),]

table(pomma_info$ci_complete)

#=== Looking for mono or poly complete CI-like prophages ===

isolates_with_ci <- summary_table %>% filter(n_complete_ci>0) 
length(unique(isolates_with_ci$isolate)) #154 lysogens with at least one phage with a complete ci
n_phages_per_isolate <- summary_table %>%  group_by(isolate) %>% tally()
n_phages_per_isolate_ci <- n_phages_per_isolate[is.element(n_phages_per_isolate$isolate,isolates_with_ci$isolate),] 
mono_poly_ci <- n_phages_per_isolate_ci %>% group_by(n) %>% tally()


