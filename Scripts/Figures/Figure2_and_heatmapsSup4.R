#Main figure 2 and heatmaps Sup fig 4

rm(list = ls())

library(ggplot2)
library(tidyr)
library(readxl)
library(edgeR)
library(magrittr)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)

colors_fig <- c("#72A86B","#A16BA8")


# ==== Fig 2A ======

gc <- read.csv("PATH_TO_GROWTH_CURVES",stringsAsFactors = F)
to_plot <- gc %>% pivot_longer(cols = starts_with("t_"), names_to = "time", names_prefix="t_", values_to="OD_600")
to_plot$time <- as.numeric(as.character(to_plot$time))

backgroud <-to_plot[which(to_plot$Isolate=="C(-)"),] 
mean_backgroud <- round(mean(backgroud$OD_600),2)

to_graph <- "AMR_A1" # AMR_A1 is the original name of CC00806 isolate

a1 <- to_plot[which(to_plot$Isolate==to_graph),]
a1$OD_600 <- a1$OD_600- mean_backgroud

a1 <- a1[which(a1$mitC!=3),]


P <- ggplot(a1,aes(time,OD_600,colour=as.factor(mitC))) +
  geom_jitter(size=0.3) + 
  scale_x_continuous(name="Time after adding MMC [h]",  breaks = seq(0,24.5,2)) +
  scale_color_manual(values=colors_fig,labels=c("Control","MMC 0.3ug/mL")) +
  ylab(bquote(OD[600])) + theme_light() +
  labs(color = "Treatment") +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 1,size = 5),legend.position = "none",
        axis.title.x = element_text(vjust = 0,size = 7),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7))


P


#=== Figure 2b===


experiment_design <- data.frame(treatment=c(rep("Control",8),rep("MMC",8)),
                                timepointh=c(rep(c(-1,0,0.5,1.25,2.25,3.25,4.25,5.5),2)),
                                score=c(rep("h",8),rep("t",8)),
                                c="c")
t <- ggplot(experiment_design, aes((timepointh), score, group = treatment, 
                                   color = treatment,shape =c ,fill=c)) +
  geom_line() + geom_point(size=3) + theme_light() +  scale_shape_manual(values = c(21)) +
  theme(axis.line = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(),panel.grid = element_blank(),
        axis.title = element_blank(),legend.position = "none",
        axis.text.y = element_blank(),axis.text.x = element_blank()) +
  scale_color_manual(values=colors_fig) +
  scale_fill_manual(values=c("white")) 

t

ggarrange(P,t,ncol = 1,hjust  = T)

#=== Figure 2C-E ===

 
load("PATH_TO_FCOUNTS_SENSE_R_OBJECT") #This was produced using Rsubread and load fc
sense <- fc
rm(fc)

mapping_file <- read.csv("PATH_TO_METADATA",
                         stringsAsFactors = F)

phages_regions <- read_xlsx("PATH_TO_ANNOTATION_PHAGES")
phages_regions$Gene <- sapply(phages_regions$Gene,function(x){strsplit(x,"AMR_A1_")[[1]][2]})
phages_regions <- phages_regions[order(phages_regions$start_position),]

singlecopy <- read.delim("PATH_TO_SINGLE_COPY_GENES",
                         stringsAsFactors = F)
rownames(singlecopy) <- singlecopy$gene


#=== parsing bam ids ===

parsing_bam_id_sense <- data.frame(library=sense$targets,
                                   sampleID=sapply(sense$targets,function(x){strsplit(x,split = "_")[[1]][2]}))

mapping_file <- merge(mapping_file,parsing_bam_id_sense,all.x = T)
mapping_file$treatment <- ifelse(mapping_file$treatment=="No mitC","control","MMC")
mapping_file$treatment <- factor(mapping_file$treatment, levels = c("control","MMC"))

group_df <- data.frame(sampleID=mapping_file$sampleID,group=paste(mapping_file$treatment,mapping_file$timepointh,sep = "_"))


#=========== Parsing counts ============

y <- featureCounts2DGEList(sense)

sample.anno <- sapply(colnames(y),function(x){strsplit(x,split = "_")[[1]][2]}) %>%
  as.data.frame() %>%
  set_colnames(c("sampleID")) # requires library(magrittr)

sample.anno <- merge(sample.anno,mapping_file,all.x = T)
sample.anno$sampleName <- sapply(sample.anno$library,function(x){strsplit(x,split = "[.]")[[1]][1]})
head(sample.anno)

group_df <- merge(sample.anno,group_df,all.x = T)
identical(group_df$sampleID,sample.anno$sampleID)

y$samples %<>%  cbind(sample.anno) %>% transform(group=factor(group_df$group,
                                                              levels = c("control_0.5", "control_1.25",
                                                                         "control_2.25", "control_3.25", "control_4.25",
                                                                         "MMC_0.5", "MMC_1.25", "MMC_2.25", "MMC_3.25",
                                                                         "MMC_4.25")))


y$samples

#=== counts ====

counts_sense <-  data.frame(library=sense$counts,strand=sense$annotation)

counts_sense$region <- "host"

for(phage in unique(phages_regions$Phage)){
  genes_phages <- NA
  genes_phages <- phages_regions$Gene[which(phages_regions$Phage==phage)]
  counts_sense$region[is.element(counts_sense$strand.GeneID,genes_phages)] <- phage
}
sum(table(counts_sense$region))
table(counts_sense$region)

sense_long <-  counts_sense %>% pivot_longer(cols = starts_with("library"),names_to = "library", names_prefix = "library." , values_to = "count" )

sense_long <- merge(sense_long,mapping_file,all.x = T)

sense_long$strand.Start1 <- as.numeric(as.character(sense_long$strand.Start))

sense_long <- sense_long[order(sense_long$strand.Start1),]

#sense_long$strand.GeneID <- factor(sense_long$strand.GeneID,levels = unique(sense_long$strand.GeneID))

sum_counts_sample <- sense_long %>% group_by(library) %>% summarise(total=sum(count))

#====== Design suggested from https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf =========

group_df <- merge(sample.anno,group_df,all.x = T)

identical(group_df$sampleID,sample.anno$sampleID)

group_df$treatment <- factor(group_df$treatment)
group_df$treatment <- relevel(group_df$treatment,ref="control")
typeof(group_df$timepointh)
group_df$timepointh <- factor(group_df$timepointh)

group_df$replicate <- factor(group_df$replicate,levels = c("3","2","1"))
design2 <- model.matrix(~replicate+treatment+timepointh+treatment:timepointh,data = group_df) # nested interaction formula
colnames(design2) <- sub("^group", "", colnames(design2))
design2
colnames(design2)


#=== Filtering to remove low counts: Comparing and deciding using filterByExpr function 

keep_edgefun <- filterByExpr(y,min.count=10,large.n=3,min.prop=1)
length(which(keep_edgefun))

ggplot(sense_long,aes(count)) + geom_histogram(binwidth = 10) + scale_x_continuous(limits = c(0,1000)) +scale_y_sqrt()
ggplot(sense_long,aes(count)) + geom_density()

before <- nrow(y$counts)
print(before)

length(which(keep_edgefun))/ before*100 #90 perc of the genes

y <- y[keep_edgefun,,keep.lib.sizes=F]
dim(y)

format(sum(keep_edgefun), big.mark = ",") # to keep genes
format(length(keep_edgefun), big.mark = ",") # all genes

sense_long$treatment2 <- as.character(sense_long$treatment)
sense_long$treatment2[!is.element(sense_long$strand.GeneID,rownames(y$counts))] <- "Filtered-out"
#ggplot(sense_long_filtered,aes(count)) + geom_density()
#ggplot(sense_long_filtered,aes(count)) + geom_histogram(binwidth = 10)  #+scale_y_sqrt()

#==== TPM ====

ggplot(sense_long,aes(count)) + geom_histogram(binwidth = 1) + scale_x_continuous(limits=c(0,100))
sense_long$rpk <- (sense_long$count)/(sense_long$strand.Length/1000)
sense_long <- sense_long %>% group_by(library) %>% mutate (tpm = rpk/(sum(rpk)/1000000))

sense_long <- sense_long %>% ungroup() %>% group_by(strand.GeneID)  %>% mutate(zscore_tpm = (log2(tpm+1) - median(log2(tpm+1)))/sd(log2(tpm+1)) ) 
# z-score based on https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-021-02936-w and https://www.biostars.org/p/107519/ 

checking_tpm <- sense_long %>% group_by(library) %>% summarise(sum(tpm))
ggplot(sense_long,aes(sampleID,asinh(tpm))) + geom_boxplot()
#checking_zscore <- sense_long %>% group_by(library) %>% summarise(sum(zscore_tpm))

phages_to_show <- sense_long[which(is.element(sense_long$region,c("LoVE","Sombra","Shia"))),]
table(interaction(phages_to_show$treatment2,phages_to_show$region))

filtered_out_phages_to_show <- phages_to_show[which(phages_to_show$treatment2=="Filtered-out"),]

exp <- ggplot(phages_to_show, aes(as.factor(timepointh),asinh(tpm),fill = treatment)) + 
  geom_violin(trim = FALSE, position = position_dodge(0.9),quantile.linewidth = 0.3,alpha=0.5) + 
  geom_boxplot(position = position_dodge(0.9),width=0.15,outlier.shape = NA, box.linewidth = 0.3,median.linewidth = 0.3,alpha=0.7)+
  geom_point(aes(colour=treatment),position = position_dodge(0.9),size=0.3,stroke = 0.3) +
  geom_hline(yintercept = asinh(100),linetype="dashed",linewidth = 0.2, color="black") +
  scale_fill_manual(values=colors_fig) +
  scale_color_manual(values = colors_fig) +
  facet_grid(vars(region)) + 
  geom_jitter(data = filtered_out_phages_to_show, colour="#4d4d4d",
            alpha=1,position = position_dodge(0.9),size=0.3) +
  scale_y_continuous(name="Asinh(TPM)") +
  scale_x_discrete(name="Time after adding MMC [h]")  +
  theme_bw()  +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 1,size = 6),legend.position = "none",
        axis.ticks.x = element_line(linewidth = 0.3),
        axis.ticks.y = element_line(linewidth = 0.3),
        axis.title.x = element_text(vjust = 0,size = 7),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.background = element_rect(fill="white"), 
        strip.text=element_text(color="black",size=7))

exp

ggarrange(P,t,exp,ncol = 2,nrow = 2,heights = c(1.5,2))


#=== Normalization for composition bias ===

y <- normLibSizes(y) #new name for calcNormFactors
head(y$samples)

#=== Exploring differences btwn libraries ===

group <- interaction(group_df$treatment,group_df$timepointh,group_df$replicate)
pch_h <- c(0,1,2,5,3,15,16,17,18,8)
colors <- rep(c("darkgreen", "red", "blue","gray","yellow"), 3)
plotMDS(y,labels = group)#col=colors[group], pch = pch_h[group])
legend("bottomright",legend = levels(group),pch = pch_h, col=colors,ncol = 3,cex = 0.5)

# mean-difference (MD) plot 

plotMD(y, column=1) # I checked all the 30 samples, all look good
abline(h=0, col="red", lty=2, lwd=2)


#=== Dispersion estimation ===

y2 <- estimateDisp(y, design2, robust=TRUE)
fit2 <- glmQLFit(y2,design2,robust = T)
head(fit2$coefficients)
summary(fit2$df.prior)
plotBCV(y2)
y2$common.dispersion  # here the common dispersion was 0.09
plotQLDisp(fit2)

#=== Differential analysis ===

colnames(design)
colnames(design2)

replicatesDiff2vs3 <- glmQLFTest(fit2, coef=2) #have to define each contrast at once
replicatesDiff1vs3 <- glmQLFTest(fit2, coef=3) #have to define each contrast at once

results_by_replicates2 <- decideTests(replicatesDiff2vs3,lfc = 2, p.value = 0.05, adjust.method = "fdr") %>% as.data.frame()
results_by_replicates1 <- decideTests(replicatesDiff1vs3,lfc = 2, p.value = 0.05, adjust.method = "fdr") %>% as.data.frame()

table(results_by_replicates2$replicate2) #there is a genuine need to adjust for the experimental replicates
table(results_by_replicates1$replicate1) #there is a genuine need to adjust for the experimental replicates

time_0.5_mitvscontrol <- glmQLFTest(fit2,coef = 4)
topTags(time_0.5_mitvscontrol)
results_replicate_t0 <- decideTests(time_0.5_mitvscontrol,lfc = 2, p.value = 0.05, adjust.method = "fdr") %>% as.data.frame()
results_replicate_t0$GeneID <- rownames(results_replicate_t0)
summary(results_replicate_t0)
table(results_replicate_t0$treatmentMMC) #All 4,889 genes not significant, so there was not difference at the base line comparison at 0.5h

time_125_mitvscontrol <- glmQLFTest(fit2,coef = 9)
topTags(time_125_mitvscontrol)
results_replicate_t1 <- decideTests(time_125_mitvscontrol,lfc = 2, p.value = 0.05, adjust.method = "fdr") %>% as.data.frame()
results_replicate_t1$GeneID <- rownames(results_replicate_t1)
summary(results_replicate_t1)
table(results_replicate_t1$treatmentMMC)  #All 4,889 genes not significant, so there was not difference at the base line comparison at 0.5h

time_225_mitvscontrol <- glmQLFTest(fit2,coef = 10)
topTags(time_225_mitvscontrol)
results_replicate_t2 <- decideTests(time_225_mitvscontrol,lfc = 2, p.value = 0.05, adjust.method = "fdr") %>% as.data.frame()
results_replicate_t2$GeneID <- rownames(results_replicate_t2)
summary(results_replicate_t2)
table(results_replicate_t2$treatmentMMC) # 42 up DEG

time_325_mitvscontrol <- glmQLFTest(fit2,coef = 11)
topTags(time_125_mitvscontrol)
results_replicate_t3 <- decideTests(time_325_mitvscontrol,lfc = 2, p.value = 0.05, adjust.method = "fdr")  %>% as.data.frame()
results_replicate_t3$GeneID <- rownames(results_replicate_t3)
summary(results_replicate_t3)
table(results_replicate_t3$treatmentMMC) # 97 up DEG

time_425_mitvscontrol <- glmQLFTest(fit2,coef = 12)
topTags(time_125_mitvscontrol)
results_replicate_t4 <- decideTests(time_425_mitvscontrol,lfc = 2, p.value = 0.05, adjust.method = "fdr")  %>% as.data.frame()
results_replicate_t4$GeneID <- rownames(results_replicate_t4)
summary(results_replicate_t4)
table(results_replicate_t4$treatmentMMC) # 103 up DEG


h1 <- as.data.frame(time_125_mitvscontrol)
h2 <- as.data.frame(time_225_mitvscontrol)
h3 <- as.data.frame(time_325_mitvscontrol)
h4 <- as.data.frame(time_425_mitvscontrol)

h1$fdr <- p.adjust(h1$PValue,method = 'fdr')
h2$fdr <- p.adjust(h2$PValue,method = 'fdr')                    
h3$fdr <- p.adjust(h3$PValue,method = 'fdr')
h4$fdr <- p.adjust(h4$PValue,method = 'fdr')

#=== Plot differential abundant per region ===

genes_df <- sense$annotation
genes_df$region <- "Host"

for(phage in unique(phages_regions$Phage)){
  genes_phages <- NA
  genes_phages <- phages_regions$Gene[which(phages_regions$Phage==phage)]
  genes_df$region[is.element(genes_df$GeneID,genes_phages)] <- phage
}
table(genes_df$region)

genes_df$filtered_out <- 0
genes_df$filtered_out[is.element(genes_df$GeneID,names(keep_edgefun)[which(!keep_edgefun)])] <- 1 
genes_df$filtered_out <- factor(genes_df$filtered_out,levels = c(0,1),labels =  c("Kept","Filtered-out")) 


genes_df <- merge(genes_df,results_replicate_t1,all.x = T)
genes_df <- merge(genes_df,results_replicate_t2,all.x = T)
genes_df <- merge(genes_df,results_replicate_t3,all.x = T)
genes_df <- merge(genes_df,results_replicate_t4,all.x = T)
genes_df$`treatmentMMC:timepointh1.25` <- factor(genes_df$`treatmentMMC:timepointh1.25`,
                                                 levels = c("Filtered-out","-1","0","1"),
                                                 labels = c("Filtered-out","Down","Not-sig","Up"))  
genes_df$`treatmentMMC:timepointh2.25` <- factor(genes_df$`treatmentMMC:timepointh2.25`,
                                                 levels = c("Filtered-out","-1","0","1"),
                                                 labels = c("Filtered-out","Down","Not-sig","Up")) 
genes_df$`treatmentMMC:timepointh3.25` <- factor(genes_df$`treatmentMMC:timepointh3.25`,
                                                 levels = c("Filtered-out","-1","0","1"),
                                                 labels = c("Filtered-out","Down","Not-sig","Up")) 
genes_df$`treatmentMMC:timepointh4.25` <- factor(genes_df$`treatmentMMC:timepointh4.25`,
                                                 levels = c("Filtered-out","-1","0","1"),
                                                 labels = c("Filtered-out","Down","Not-sig","Up"))    

genes_df_long <- genes_df %>% select(GeneID,region,filtered_out,`treatmentMMC:timepointh1.25`,
                                     `treatmentMMC:timepointh2.25`,`treatmentMMC:timepointh3.25`,
                                     `treatmentMMC:timepointh4.25`) %>% pivot_longer(!c(GeneID,region,filtered_out))  

genes_df_long$time_point <- sapply(genes_df_long$name,function(x){strsplit(x,split = "inth")[[1]][2]})

identical(unique(genes_df_long$GeneID[which(is.na(genes_df_long$value))]),
          unique(genes_df_long$GeneID[which(genes_df_long$filtered_out=="Filtered-out")]))

genes_df_long$value[which(is.na(genes_df_long$value))] <- "Filtered-out"

genes_not_significant_differences <- genes_df[which(genes_df$`treatmentMMC:timepointh1.25`=="Not-sig"&
                                                      genes_df$`treatmentMMC:timepointh2.25`=="Not-sig"&
                                                      genes_df$`treatmentMMC:timepointh3.25`=="Not-sig"&
                                                      genes_df$`treatmentMMC:timepointh4.25`=="Not-sig"),]

n_genes_per_region <- genes_df %>% group_by(region) %>% summarise(n_genes_per_region=n())

sum_per_region <- genes_df_long %>% group_by(region,time_point,value) %>% summarise(n_genes=n())
sum_per_region <- merge(sum_per_region,n_genes_per_region,all.x = T)

sum_per_region$region <- factor(sum_per_region$region,levels = c("Host","LoVE","Shia","Sombra","Hanky","Wilby_1","Wilby_2"))



g <- ggplot(sum_per_region[which(sum_per_region$region!="Host"),],aes(time_point,n_genes/n_genes_per_region*100,fill=value)) + 
  scale_fill_manual(values=brewer.pal(3,"Greys"),
                    labels=c("Filtered-out","Not-sig","Up-regulated")) +
  geom_bar(stat = "identity") + facet_grid(region~.) + theme_light() +
  scale_y_continuous(name = "Percentage of genes per DA comparison category",limits=c(0,100),breaks = seq(0,100,50)) +
  scale_x_discrete(name="Time of MMC vs control comparision [h]") +
  #labs(fill = "") +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 1,size = 6),legend.position = "none",
        axis.title.x = element_text(vjust = 0,size = 7),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.background = element_rect(fill="white"), 
        strip.text=element_text(color="black",size=7))
g


left <- ggarrange(P,t,g,ncol = 1, nrow = 3,heights = c(1.2,1,2))
left

# === Setting colors heatmap ====

max(asinh(sense_long$tpm))
dens <- density(asinh(sense_long$tpm))
df_dense <- data.frame(x=dens$x, y=dens$y)

probs <- seq(0,1,by=0.1)
quantiles <- quantile(df_dense$x, prob=probs)

quantile(df_dense$x,0.3) #3.2424356 starts the second peak in the kernel density estimation

sinh(3.2424356) # 12.8 TPM have a second peak 
asinh(12.8)

3.012479 # the peak really starts in this point(just after the lowest value )
sinh(3.012479) # 10.14 TPM have a second peak 
asinh(10.14)


df_dense$quant <- factor(findInterval(df_dense$x,quantiles))

brewer.pal(9,"Greys")
cols <- colorRampPalette(brewer.pal(9,"OrRd"))(length(quantiles))
names(cols) <- levels(df_dense$quant)
length(cols)

ggplot(df_dense, aes(x,y)) + geom_line() + geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) +
  scale_x_continuous(breaks=(quantiles)) + scale_fill_manual(values=cols)

ggplot(df_dense, aes(x,y)) + geom_line() + geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) +
  scale_x_continuous() + scale_fill_manual(values=cols)

quantiles


# === heatmap ===

phage <- "LoVE"

str(sense_long)
sense_long$timepointh <- factor(as.character(sense_long$timepointh),levels = as.character(c(0.5,1.25,2.25,3.25,4.25)))
sense_long <- sense_long[order(sense_long$treatment,sense_long$timepointh),]
unique(sense_long$sampleID)
sense_long

give_pheatmap <- function(phage){
  
  z_p_heatmap <- sense_long %>% filter(region==phage)  %>% pivot_wider(names_from = library,values_from = tpm, id_cols = strand.GeneID)
  z_p_heat <- as.matrix(z_p_heatmap[,2:ncol(z_p_heatmap)])
  rownames(z_p_heat) <- z_p_heatmap$strand.GeneID
  
  #== treatment
  
  mp_zp_heat <- sense_long %>% ungroup %>% group_by(treatment,timepointh,replicate,library) %>% summarise(n=n())
  z_p_heatmap <- z_p_heatmap[,c("strand.GeneID",mp_zp_heat$library)]
  identical(mp_zp_heat$library,colnames(z_p_heatmap)[-1])
  
  treatment_df <- data.frame(treatment=mp_zp_heat$treatment)
  row.names(treatment_df) <- colnames(z_p_heatmap)[-1] 
  
  #== single copy genes
  
  sc_love <- singlecopy[z_p_heatmap$strand.GeneID,]
  identical(sc_love$gene,(z_p_heatmap$strand.GeneID))
  sc_df <- data.frame(Single_copy=sc_love$single_copy)
  rownames(sc_df) <- sc_love$gene
  
  col_treatments <- colors_fig
  names(col_treatments) <- c("control","MMC")
  
  
  labels_love <- phages_regions[which(phages_regions$Phage==phage),] %>% as.data.frame()
  
  if(phage=="LoVE"){
    labels_love$Manual[which(labels_love$Gene=="contig00001_2120")] <- "cysH"
  }
  
  if(phage=="Sombra"){
    labels_love <- labels_love[-grep("rRNA",labels_love$Gene),]
  }
  identical(labels_love$Gene,z_p_heatmap$strand.GeneID)
  
  labels_love$Manual[which(labels_love$Manual=="Hypothetical")]<- ""
  
  col_sc <- c("#FCF6F5","#61bdcd")#"#2BAE66") #https://designwizard.com/blog/colour-combination/
  names(col_sc) <- c("Yes","No")
  colnames(sc_df)[1] <- "sc"
  mycolors2 <- list(treatment=col_treatments,sc=col_sc)
  
  if(phage=="LoVE"){
    #Adding DA 
    da <- genes_df %>% select(`treatmentMMC:timepointh4.25`,`treatmentMMC:timepointh3.25`,
                              `treatmentMMC:timepointh2.25`,`treatmentMMC:timepointh1.25`) %>% as.data.frame()
    da[is.na(da)] <- "Filtered-out"
    rownames(da) <- genes_df$GeneID 
    da_love <- da[z_p_heatmap$strand.GeneID,]
    identical(rownames(da_love),(z_p_heatmap$strand.GeneID))
    
    sc_df <- cbind(sc_df,da_love)
    colnames(sc_df)<- c("sc","t4","t3","t2","t1")
    
    identical(rownames(sc_df),labels_love$Gene)
    #sc_df$gene_cat <- labels_love$Manual_category
    
    brewer.pal(3,"Greys")
    
    col_t1 <- c("#F0F0F0","#BDBDBD")#,"#a8b197")
    names(col_t1) <- c("Filtered-out","Not-sig")
    
    col_t2 <- c("#F0F0F0", "#BDBDBD" ,"#636363")#,"#a8b197")
    names(col_t2) <- c("Filtered-out","Not-sig","Up")
    mycolors2 <- list(treatment=col_treatments,sc=col_sc,t1=col_t1,
                      t2=col_t2,t3=col_t2,t4=col_t2)
    
  }
  
  
  #use the asinh to normalise
  h <- pheatmap(asinh(z_p_heat),cluster_cols = F,cluster_rows = F,
                breaks = quantiles, color = cols,
                annotation_col = treatment_df,
                #annotation_col_legend_width=1,
                annotation_row = sc_df,
                labels_row = labels_love$Manual,
                show_colnames = F, show_rownames = T,
                fontsize_row = 5,
                #fontsize_col = 5,
                fontsize=5,
                #fontsize = 5,
                cellwidth = 3, cellheight = 3,
                gaps_col = c(3,6,9,12,15,15,18,21,24,27),
                annotation_colors = mycolors2,
                trace="none",
                legend = T, 
                annotation_legend = F) 
  
  return(h)
  
}

h <- give_pheatmap(phage)

left <- ggarrange(P,t,g,ncol = 1, nrow = 3,heights = c(0.9,1,2))
left

ggsave(left, width = 6, height = 16 , units = 'cm', dpi = 320,
     filename = "PATH_OUTPUT_LEFT_SIDE_FIG2")


ggsave(h, width = 12 , height = 16 , units = 'cm', dpi = 320,
       filename = "PATH_OUTPUT_LOVE_HEATMAP")


#=== Sup Fig 4: Other Phages "Hanky    host    LoVE    Shia  Sombra Wilby_1 Wilby_2 "====

w1 <- give_pheatmap("Wilby_1")
w2 <- give_pheatmap("Wilby_2")

sombra <- give_pheatmap("Sombra")
hanky <- give_pheatmap("Hanky")
shia <- give_pheatmap("Shia")


ggsave(w1, width =  9 , height = 6 , units = 'cm', dpi = 320,
       filename = "PATH_OUTPUT_W1_HEATMAP")

ggsave(w2, width =  9 , height = 7, units = 'cm', dpi = 320,
       filename = "PATH_OUTPUT_W2_HEATMAP")

ggsave(shia, width =  9 , height = 9 , units = 'cm', dpi = 320,
       filename = "PATH_OUTPUT_SHIA_HEATMAP")

ggsave(sombra, width =  9 , height = 9 , units = 'cm', dpi = 320,
       filename = "PATH_OUTPUT_SOMBRA_HEATMAP")

ggsave(hanky, width =  9 , height = 6 , units = 'cm', dpi = 320,
       filename = "PATH_OUTPUT_HANKY_HEATMAP")



