rm(list=ls())

library(rhmmer) #https://github.com/arendsee/rhmmer
library(ggplot2)
library(dplyr)

#========= Loading files ======

jackhmmer_domains <- rhmmer::read_domtblout("PATH_TO_DOMAIN_TABLE")


#=== Getting coverage ===

jackhmmer_domains$hmm_coverage <- jackhmmer_domains$hmm_to - jackhmmer_domains$hmm_from + 1 

ggplot(jackhmmer_domains,aes(hmm_coverage)) + geom_histogram()

jackhmmer_domains$align_protein_coverage <- jackhmmer_domains$env_to - jackhmmer_domains$env_from +1 
jackhmmer_domains$protein_coverage <- jackhmmer_domains$align_protein_coverage/jackhmmer_domains$domain_len *100

hist(jackhmmer_domains$protein_coverage)

quantile(jackhmmer_domains$domain_ievalue)
quantile(jackhmmer_domains$domain_score)

#=== Filtering out hits with domain bitscore>=30 and protein_coverage (this is the phage protein) >=70 ===

jackhmmer_domains <- jackhmmer_domains[which( jackhmmer_domains$domain_score>=30 &jackhmmer_domains$protein_coverage>=70 ),]

hist(jackhmmer_domains$protein_coverage)

#=== Checking scores ===

quantile(jackhmmer_domains$domain_score)
quantile(jackhmmer_domains$domain_ievalue)
quantile(jackhmmer_domains$sequence_evalue)
quantile(jackhmmer_domains$sequence_score)

#===  Are domains repeated? No in this case ===

n_times_protein_per_iteration <- jackhmmer_domains %>% group_by(domain_name) %>% summarise(n=n())

#==== Counting subgroups:Calculating number of proteins in each quadrant ===

jackhmmer_domains$subgroup <- ifelse(jackhmmer_domains$hmm_coverage>=237*0.7,"Cov>=70","Cov<70")

table(jackhmmer_domains$subgroup) #1,666 cov>=70 and 442 cov<70

#==== Adding region which match to ===

#Carboxy terminal domain

jackhmmer_domains$domain <- NA

# Domain boundary positions (adjust to your HMM)
pos_first_helix     <- 33    # Position of first helix (e.g. Gln Q34)
pos_active_site     <- 193   # Position of active site (e.g. Lys193)
pos_end_hinge       <- 132   # End of hinge region
ntd_truncated_lim   <- 92    # End of the NTD

for(i in 1:nrow(jackhmmer_domains)){
  
  if(jackhmmer_domains$hmm_from[i] <= pos_first_helix){ #position gln Q 34 first helix
    
    if(jackhmmer_domains$hmm_to[i] > pos_active_site) { #position lys193 active site
      
      jackhmmer_domains$domain[i] <- "whole_protein"
      
    }else if(jackhmmer_domains$hmm_to[i] <= pos_end_hinge) { # end hinge
      
      jackhmmer_domains$domain[i] <- "NTD"
      
    } else{
      
      jackhmmer_domains$domain[i] <- "NTD-and-CTD_truncated"
      
    }
    
  }else if (jackhmmer_domains$hmm_from[i] <= ntd_truncated_lim){
    
    if(jackhmmer_domains$hmm_to[i] > pos_active_site) { #position lys193 active site
      
      jackhmmer_domains$domain[i] <- "NTD_truncated_and_CDT"
      
    }else if(jackhmmer_domains$hmm_to[i] <= pos_end_hinge) { # end hinge
      
      jackhmmer_domains$domain[i] <- "NTD_truncated"
      
    } else{
      
      jackhmmer_domains$domain[i] <- "NTD_truncated-and-CTD_truncated"
      
    }
  }else{
    jackhmmer_domains$domain[i] <- "CTD_truncated"
  }

  
  if(is.na(jackhmmer_domains$domain[i])){
    jackhmmer_domains$domain[i] <- "check"
  }
  
  
}

table(jackhmmer_domains$domain,useNA = "ifany")
ggplot(jackhmmer_domains,aes(domain,domain_len)) + geom_boxplot()

table(jackhmmer_domains$domain[which(startsWith(jackhmmer_domains$domain_name,"bacteroidales"))])

417-199 # 218 from non-Bacteroidales isolates

write.table(jackhmmer_domains,"PATH_TO_PARSED_JACKHMMER_RESULTS",
            quote = F, sep = "\t",row.names = F)
