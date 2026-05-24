#Colors tree phylo CI complete


rm(list=ls())
library(tidyverse)


#=== Loading files ===

phages_no_bacteroidales <- read.delim("PATH_TO_SUMMARY_NO_BACTEROIDALES_PROPHAGES")

phages_bacteroidales <- read.delim("PATH_TO_SUMMARY_BACTEROIDALES_PROPHAGES")

tree_proteins <- read.delim("PATH_TO_ALL_PROTEINS_IN_THE_TREE",
                            header = F)

#=== Getting hosts' taxa info ====

tree_proteins$PhageID <- gsub("_[0-9]*$","",tree_proteins$V1)
tree_proteins$Phylum <- NA

for(i in 1:nrow(tree_proteins)){
  if(startsWith(tree_proteins$V1[i],"bacteroidales")){
    tree_proteins$Phylum[i] <- "Bacteroidota"
  } else if (tree_proteins$V1[i]=="sp|P03034|RPC1_LAMBD"){
    tree_proteins$Phylum[i] <- "Pseudomonadota"
  } else {
    tree_proteins$Phylum[i] <- phages_no_bacteroidales$phylum[phages_no_bacteroidales$vs2_2nd_round_name==tree_proteins$PhageID[i]]
  }
}

#=== Getting leave ID ===

tree_proteins$leave <- gsub(":","_",tree_proteins$V1)
tree_proteins$leave <- gsub("\\|","_",tree_proteins$leave)

#=== Defining colors same as Sof's ms ==

colors <- data.frame(Phylum=unique(tree_proteins$Phylum),colour=NA)
colors$colour[which(colors$Phylum=="Bacillota")] <- "#b7cae0"##0b5394ff"
colors$colour[which(colors$Phylum=="Pseudomonadota")] <- "#fac2b9" #"#cc0000ff"
colors$colour[which(colors$Phylum=="Bacteroidota")] <- "#ecdec1"#"#a3b59d" #"#274e13ff"
colors$colour[which(colors$Phylum=="Actinomycetota")] <- "#bac7b5"##bf9000ff"
colors$colour[which(colors$Phylum=="Fusobacteriota")] <- "#9e9e9e"#000000ff"


tree_proteins <- merge(tree_proteins,colors)


to_print_range <- tree_proteins %>% select(leave,colour,Phylum) %>% 
  mutate(range="range", .before = colour) 

to_print_range$colour[grep("sp",to_print_range$leave)] <- "#cc0000ff"
to_print_range$Phylum[grep("sp",to_print_range$leave)] <- "Pseudomonadota (Lambda CI)"


write.csv(to_print_range,file = "PATH_TO_COLOUR_RANGE_PHYLA_OUTPUT",
          row.names = F, quote = F)

#cat ../../Figures/Figure1/info_Itol_tree/template_range.txt colour_strip_phyla_hos_LINSI_conservative.csv > itol_range_host_phyla_SupFig1.txt

#=== Complete CI genomic location ===


complete_ci <- phages_bacteroidales[which(phages_bacteroidales$n_complete_ci>0),]
complete_ci$orf_ci <- NA
complete_ci$per_position <- NA
complete_ci$gene_count <- as.numeric(as.character(complete_ci$gene_count))

for(i in 1:nrow(complete_ci)){
  
  if(grepl(";",complete_ci$complete_ci[i])){
    proteins <- unlist(strsplit(complete_ci$complete_ci[i],split = ";"))
    protein_name_array <- sapply(proteins, function(x){strsplit(x,split = "_")})
    proteins_orf <- c()
    for(j in 1:length(protein_name_array)){
      proteins_orf <- c(proteins_orf,protein_name_array[[j]][length(protein_name_array[[j]])])
    }
    
    complete_ci$orf_ci[i] <- paste(proteins_orf,collapse = ";")
    complete_ci$per_position[i] <- max(as.numeric(as.character(proteins_orf))/complete_ci$gene_count[i])
    
  }else{
    protein_name_array <- strsplit(complete_ci$complete_ci[i],split = "_")
    complete_ci$orf_ci[i] <- protein_name_array[[1]][length(protein_name_array[[1]])]
    complete_ci$per_position[i] <- as.numeric(as.character(complete_ci$orf_ci[i]))/complete_ci$gene_count[i]
  }
}



p <- ggplot(complete_ci,aes(per_position)) + geom_histogram(bins = 15) +
  theme_light() +  xlab("Relative position in the prophage genome") +
  ylab("N phages with complete-CI in this position") +
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 1,size = 5),
        axis.title.x = element_text(vjust = 0,size = 7),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7)) 


ggsave(p, width = 7.5, height = 7 , units = 'cm', dpi = 320,
       filename ="PATH_TO_SUP_FIG_1B")


ggplot(complete_ci,aes(per_position)) + geom_density()


#== Complete CI

ci_no_bacte <- phages_no_bacteroidales %>% filter(n_complete_ci>0)
table(ci_no_bacte$phylum) 

ci_bacte <- phages_bacteroidales %>% filter(n_complete_ci>0)
table(ci_bacte$phylum) 


ntd_no_bacte <- phages_no_bacteroidales %>% filter(n_ntd_ci>0)
table(ntd_no_bacte$phylum) 

ntd_bacte <- phages_bacteroidales %>% filter(n_ntd_ci>0)
table(ntd_bacte$phylum) 



