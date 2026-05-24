#Comparing Wilbys

#=== Organizing the environment ===

rm(list=ls())
library(readxl)
library(gggenes)

#==== Loading the files =====

a1_genome <- read.delim("PATH_TO_CC00806_ANNOTATIONS")

phages_regions <- read_xlsx("PATH_TO_PHAGE_ANNOTATIONS")
phages_regions$Gene <- sapply(phages_regions$Gene,function(x){strsplit(x,"AMR_A1_")[[1]][2]})
phages_regions <- phages_regions[order(phages_regions$start_position),]

singlecopy <- read.delim("PATH_TO_SINGLE_COPY_PROTEINS",
                         stringsAsFactors = F) #it came of a cdhit at 70% global 
colnames(singlecopy)[1] <- "Gene"

#=== Taking Wilby ===

wilby_regions <- phages_regions[is.element(phages_regions$Phage,c("Wilby_1","Wilby_2")),]

wilby_regions <- merge(wilby_regions,singlecopy,all.x = T) 

wilby_regions <- wilby_regions %>% select(Phage,cluster,Gene,start_position,end_position,strandedness,Manual_category,Manual)

wilby_regions <- wilby_regions %>% group_by(cluster) %>% mutate(n=n())

wilby_regions <- wilby_regions[which(wilby_regions$Gene!="contig00001_2443"),]

wilby_regions <- wilby_regions %>% group_by(Phage) %>% mutate(min_loc=min(start_position),end_loc=max(end_position))

wilby_regions$relative_start_position_kb <- (wilby_regions$start_position - wilby_regions$min_loc)
wilby_regions$relative_end_position_kb <- (wilby_regions$end_position - wilby_regions$min_loc)

wilby_regions$orientation <- ifelse(wilby_regions$strandedness==1,T,F)

#There is no alignments between 24,614bp and 25,870 bp in Wilby_1 (query)  and between 25,656 bp to 27,141 bp in Wilby_2

df_no_align <- data.frame()

alignments <- data.frame(x=c(693,24613,25655,1728,693,
                             24661,25057,26180,25783,24661,
                             25871,43128,44398,27142,25871,
                             652,1178,1726,1203,652,
                             40975,41086,1120,1009,40967),
                         y=rep(c(rep("Wilby_1",2),rep("Wilby_2",2),"Wilby_1"),5),
                         group=c(rep("a",5),rep("b",5),rep("c",5),rep("d",5),rep("e",5)))

genes_to_color <- c("contig00001_1064","contig00001_1065","contig00001_1100","contig00001_1101","contig00001_1111",
                    "contig00001_2444","contig00001_2445","contig00001_2446","contig00001_2485","contig00001_2486","contig00001_2498")

wilby_regions$Manual_category[!is.element(wilby_regions$Gene,genes_to_color)] <- "other"
wilby_regions$Manual_category <- factor((wilby_regions$Manual_category),
       levels=c("Hypothetical","CI-like protein","moron, auxiliary metabolic gene and host takeover", "integration and excision", "other"))

wilby_regions$Manual[!is.element(wilby_regions$Gene,genes_to_color)] <- NA

df_colors <- data.frame(categories=factor(unique(wilby_regions$Manual_category),
                                          levels=c("Hypothetical","CI-like protein",
                                                   "moron, auxiliary metabolic gene and host takeover",
                                                   "integration and excision",
                                                   "other")), colors = NA)

df_colors <- df_colors[order(df_colors$categories),]

df_colors$colors <- c("#bdbdbd","#e31a1c","#ff7f00","#fb9a99","white")




p <- ggplot(wilby_regions, aes(xmin = relative_start_position_kb, xmax = relative_end_position_kb, y=Phage,forward=orientation,fill=Manual_category)) +
  geom_text_repel(data=wilby_regions %>% mutate(start2 = (relative_start_position_kb+relative_end_position_kb)/2), 
                  aes(x=start2, label = Manual),  size = 5 / .pt,
                  angle= 80,nudge_y = 0.1,
                  segment.linetype = "dashed",
                  segment.size = 0.1) + coord_cartesian(clip = "off")  +
  scale_fill_manual(values = df_colors$colors) +
  geom_polygon(data = alignments[which(alignments$group!="e"),],aes(x=x,y=y,group=(group)),inherit.aes = F,alpha=0.8)+
  geom_gene_arrow(arrowhead_width = unit(2, "mm")) +
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
                                  #axis.ticks.x = element_blank(),
                                  #plot.margin = unit(c(1,1,1,1),"lines"))+
                                  panel.grid.major = element_blank(), 
                                  panel.grid.minor = element_blank()) +
  xlab("Prophage genomic location [bp]") + ylab("")

  
p


ggsave(p, width = 18, height = 10, units = 'cm', dpi = 320,
       filename = "PATH_TO_SUP_FIG_3_WILBIES")



#=== surrounding areas ===

genes_to_fig  <- c("AMR_A1_contig00001_1061","AMR_A1_contig00001_1062","AMR_A1_contig00001_1063",
           "AMR_A1_contig00001_1112","AMR_A1_contig00001_1113","AMR_A1_contig00001_1114",
           "AMR_A1_contig00001_2441","AMR_A1_contig00001_2442","AMR_A1_contig00001_2443",
           "AMR_A1_contig00001_2499","AMR_A1_contig00001_2500","AMR_A1_contig00001_2501")

fig_surr <- a1_genome[is.element(a1_genome$X,genes_to_fig),]

fig_surr$Phage <- "Wilby_1"

fig_surr$Phage[is.element(fig_surr$X, c("AMR_A1_contig00001_2441","AMR_A1_contig00001_2442","AMR_A1_contig00001_2443",
                                           "AMR_A1_contig00001_2499","AMR_A1_contig00001_2500","AMR_A1_contig00001_2501"))] <- "Wilby_2"

fig_surr$Gene <- sapply(fig_surr$X,function(x){strsplit(x,split = "A1_")[[1]][2]})

fig_surr <- merge(fig_surr,singlecopy,all.x = T)



fig_surr <- fig_surr %>% select(Phage,Gene,start_position,end_position,cluster,kegg_hit,pfam_hits)

wilbys_ends <- wilby_regions %>% group_by(Phage) %>% summarise(start_position=min(start_position),end_position=max(end_position))

wilbys_ends <- data.frame(Phage=wilbys_ends$Phage,Gene=wilbys_ends$Phage,start_position=wilbys_ends$start_position,
                          end_position=wilbys_ends$end_position,cluster=NA,kegg_hit=NA,pfam_hits=NA)

round((wilbys_ends$end_position- wilbys_ends$start_position+1)/1000,1)

wilbys_dS_tRNA <- data.frame(Phage=c("Wilby_1","Wilby_2"),Gene=c("tRNA_1","tRNA_2"),start_position=c(1481488,3137814),
                          end_position=c(1481416,3137741),cluster=NA,kegg_hit=NA,pfam_hits=NA)

fig_surr <- rbind(fig_surr,wilbys_ends,wilbys_dS_tRNA)

p <- ggplot(fig_surr[which(fig_surr$Phage=="Wilby_1"),], aes(xmin = start_position, xmax = end_position, y=Phage, fill = Gene)) +
  geom_gene_arrow(arrowhead_width = unit(2, "mm")) + facet_grid(Gene~.,scales  = "free")


p

# trna downstream phages 1)/  1481488 1481416 Trp     CCA     72.9 and 2)    69      3137814 3137741 Arg     TCT     72.1
