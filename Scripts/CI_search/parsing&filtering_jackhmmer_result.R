rm(list=ls())
library(rhmmer)  # https://github.com/arendsee/rhmmer
library(ggplot2)
library(dplyr)

# ========= PATHS =========

#input_file  <- "<PATH_TO_JACKHMMER_DOMTBLOUT_FILE>"   # e.g. all_phages_proteome_dom_7_iter.txt
#output_file <- "<PATH_TO_OUTPUT_TSV_FILE>"             # e.g. parsed_jackhmmer_all_phages.tsv


# ========= PARAMETERS =========
min_domain_score    <- 30    # Minimum domain bitscore
min_protein_coverage <- 70   # Minimum phage protein coverage (%)
hmm_total_length    <- 237   # Total HMM length (used for HMM coverage subgrouping)
hmm_subgroup_cutoff <- 0.7   # Fraction of HMM length for Cov>=70 subgroup

# Domain boundary positions (adjust to your HMM)
pos_first_helix     <- 33    # Position of first helix (e.g. Gln Q34)
pos_active_site     <- 193   # Position of active site (e.g. Lys193)
pos_end_hinge       <- 132   # End of hinge region

# ========= Loading files =========
jackhmmer_domains <- rhmmer::read_domtblout(input_file)

# ========= Getting coverage =========
jackhmmer_domains$hmm_coverage <- jackhmmer_domains$hmm_to - jackhmmer_domains$hmm_from + 1
ggplot(jackhmmer_domains, aes(hmm_coverage)) + geom_histogram()

jackhmmer_domains$align_protein_coverage <- jackhmmer_domains$env_to - jackhmmer_domains$env_from + 1
jackhmmer_domains$protein_coverage <- jackhmmer_domains$align_protein_coverage / jackhmmer_domains$domain_len * 100
hist(jackhmmer_domains$protein_coverage)

quantile(jackhmmer_domains$domain_ievalue)
quantile(jackhmmer_domains$domain_score)

# ========= Filtering: domain bitscore >= threshold and protein coverage >= threshold =========

jackhmmer_domains <- jackhmmer_domains[which(
  jackhmmer_domains$domain_score    >= min_domain_score &
    jackhmmer_domains$protein_coverage >= min_protein_coverage
), ]
hist(jackhmmer_domains$protein_coverage)

# ========= Checking scores =========
quantile(jackhmmer_domains$domain_score)
quantile(jackhmmer_domains$domain_ievalue)
quantile(jackhmmer_domains$sequence_evalue)
quantile(jackhmmer_domains$sequence_score)

# ========= Are domains repeated? =========
n_times_protein_per_iteration <- jackhmmer_domains %>%
  group_by(domain_name) %>%
  summarise(n = n())

# ========= Counting subgroups by HMM coverage =========
jackhmmer_domains$subgroup <- ifelse(
  jackhmmer_domains$hmm_coverage >= hmm_total_length * hmm_subgroup_cutoff,
  "Cov>=70", "Cov<70"
)
table(jackhmmer_domains$subgroup)

# ========= Assigning domain regions =========
jackhmmer_domains$domain <- NA

for (i in 1:nrow(jackhmmer_domains)) {
  
  if (jackhmmer_domains$hmm_from[i] <= pos_first_helix) {
    
    if (jackhmmer_domains$hmm_to[i] > pos_active_site) {
      jackhmmer_domains$domain[i] <- "whole_protein"
      
    } else if (jackhmmer_domains$hmm_to[i] <= pos_end_hinge) {
      jackhmmer_domains$domain[i] <- "NTD"
      
    } else {
      jackhmmer_domains$domain[i] <- "NTD-and-CTD_truncated"
    }
    
  } else if (jackhmmer_domains$hmm_from[i] <= pos_end_hinge) {
    
    if (jackhmmer_domains$hmm_to[i] > pos_active_site) {
      jackhmmer_domains$domain[i] <- "NTD_truncated_and_CDT"
      
    } else if (jackhmmer_domains$hmm_to[i] <= pos_end_hinge) {
      jackhmmer_domains$domain[i] <- "NTD_truncated"
      
    } else {
      jackhmmer_domains$domain[i] <- "NTD_truncated-and-CTD_truncated"
    }
    
  } else {
    jackhmmer_domains$domain[i] <- "CTD_truncated"
  }
  
  if (is.na(jackhmmer_domains$domain[i])) {
    jackhmmer_domains$domain[i] <- "check"
  }
}

table(jackhmmer_domains$domain, useNA = "ifany")
ggplot(jackhmmer_domains, aes(domain, domain_len)) + geom_boxplot()

# ========= Optional: filter by taxonomy prefix =========
taxonomy_prefix <- "<TAXONOMY_PREFIX>"   # e.g. "bacteroidales"
table(jackhmmer_domains$domain[which(startsWith(jackhmmer_domains$domain_name, taxonomy_prefix))])

# ========= Export =========
write.table(jackhmmer_domains, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
