library(ggplot2)
library(biomaRt)
library(dplyr)

# Read ASE data from SplAdder database
ASE_Spladder<- read.table("/Users/aychaserradj/Desktop/ZHAW/Praktikum/data/SplAdder_PSI_COAD_and_READ.tsv", header = TRUE, sep = "\t")
ASE_Spladder$GeneID <- sub("\\..*", "", ASE_Spladder$GeneID)
ASE_SpliceSeq<- read.table("/Users/aychaserradj/Desktop/ZHAW/Praktikum/data/PSI_COAD_and_READ.tsv", header = TRUE, sep = "\t")

# rename ASE types to common names in SplAdder and SpliceSeq data base
# Define a named vector with the replacements
replacements_spliceseq <- c("AA" = "Alternate acceptor",
                            "AD" = "Alternate donor",
                            "ES" = "Exon skipping",
                            "ME" = "Mutually exclusive exons",
                            "RI" = "Retained Intron",
                            "AP" = " Alternate promotor",
                            "AT" = "Alternate terminator")

replacements_spladder <- c("A3" = "Alternate acceptor", 
                           "A5" = "Alternate donor",
                           "ES" = "Exon skipping",
                           "ME" = "Mutually exclusive exons",
                           "IR" = "Retained Intron")

ASE_Spladder$SpliceType <- replacements_spladder[ASE_Spladder$SpliceType]
ASE_SpliceSeq$splice_type <- replacements_spliceseq[ASE_SpliceSeq$splice_type]

# Get column names that end with the specific substring
# for SplAdder data base
columns_with_substring_normal <- grep("Event$|Type$|Symbol$|GeneID$|.11A.", names(ASE_Spladder), value = TRUE)
columns_with_substring_tumor <- grep("Event$|Type$|Symbol$|GeneID$|.01A.", names(ASE_Spladder), value = TRUE)

# Create a subset with tumor and normal samples
normal_samples <- subset(ASE_Spladder, select = columns_with_substring_normal)
tumor_samples <- subset(ASE_Spladder, select = columns_with_substring_tumor)

# Select normal samples of SpliceSeq data
culumns_spliceseq_normal <- grep("symbol$|as_id$|splice_type$|_Norm",  names(ASE_SpliceSeq), value = TRUE)
normal_samples_spliceseq <- subset(ASE_SpliceSeq, select = culumns_spliceseq_normal)
# Pattern for tumor samples
#  Select tumor samples of SpliceSeq data
pattern <- "^[^_]*_[^_]*_[^_]*$|symbol$|as_id$|splice_type$"
selected_columns <- grep(pattern, names(ASE_SpliceSeq), value = TRUE)

# Subset the data frame with selected columns
tumor_samples_spliceseq <- ASE_SpliceSeq[selected_columns]
tumor_samples_spliceseq <- tumor_samples_spliceseq[, -3]


# source: https://stat.ethz.ch/pipermail/r-help/2008-February/154167
perform_t_test <- function(t, n) {
  obj <- try(t.test(t, n), silent = TRUE)
  if (is(obj, "try-error")) return (NA) else return(obj$p.value)
}


# Method to perform t-test
run_t_test <- function(normal_samples, tumor_samples, db){
  results_p_values <- vector(mode = "numeric", length = 0)
  ASE_gene_names <- character(0)
  ASE_type <- character(0)
  ASE_names <- character(0)
  for (i in 1:nrow(normal_samples)){
    tumor <- as.numeric(tumor_samples[i,])
    normal <- as.numeric(normal_samples[i,])
    # for spliceseq only remove first two one, SpliceAdder remove first tree one
    if (db == "spliceseq"){
      tumor_t <- tumor[-c(1,2,3)]
      normal_t <- normal[-c(1,2,3)]
    } else if (db == "spladder"){
      tumor_t <- tumor[-c(1,2,3,4)]
      normal_t <- normal[-c(1,2,3,4)]
    }

    if (db == "spladder"){
      ASE_type <- c(ASE_type, paste(normal_samples[i,][4]))
      ASE_names <- c(ASE_names, paste(normal_samples[i,][1]))
      ASE_gene_names <- c(ASE_gene_names, paste(normal_samples[i,][2]))
    }else{
      ASE_type <- c(ASE_type, paste(normal_samples[i,][3]))
      ASE_names <- c(ASE_names, paste(normal_samples[i,][2]))
      ASE_gene_names <- c(ASE_gene_names, paste(normal_samples[i,][1]))
    }
    t_test_p_value <- perform_t_test(tumor_t, normal_t)
    results_p_values <- c(results_p_values, t_test_p_value)
  }
  # Create dictionary of ASE names and p_values from t-test
  results_p_values <- p.adjust(results_p_values, method = "BH")
  
  indices_to_remove <- c()
  
  # Iterate over the indices of results_p_values
  indices_to_remove <- which(is.na(results_p_values) | is.nan(results_p_values) | results_p_values >= 0.05)
  
  # Remove elements from all vectors simultaneously
  results_p_values <- results_p_values[-indices_to_remove]
  ASE_type <- ASE_type[-indices_to_remove]
  ASE_names <- ASE_names[-indices_to_remove]
  ASE_gene_names <- ASE_gene_names[-indices_to_remove]
  
  ASE_and_p_values <- setNames(ASE_gene_names, results_p_values)
  return(list(ASE_and_p_values, ASE_type, ASE_gene_names))
}


SpliceSeq_ASE_types <- as.data.frame(ASE_SpliceSeq_p_values[2])
SplAdder_ASE_types <- as.data.frame(ASE_SplAdder_p_values[2])
names(SpliceSeq_ASE_types)[names(SpliceSeq_ASE_types) == "c..Retained.Intron....Retained.Intron....Alternate.donor....Exon.skipping..."] <- "ase_type"

SpliceSeq_ASE_types_percentage <- prop.table(table(SpliceSeq_ASE_types$ase_type)) * 100
SplAdder_ASE_types_percentage <- prop.table(table(SplAdder_ASE_types$ase_type)) * 100

SplAdder_ASE_types_percentage <- as.data.frame(SplAdder_ASE_types_percentage)
SpliceSeq_ASE_types_percentage <- as.data.frame(SpliceSeq_ASE_types_percentage)

SplAdder_ASE_types_percentage$group <- "SplAdder"
SpliceSeq_ASE_types_percentage$group <- "SpliceSeq"

ASE_types_combined <- bind_rows(SplAdder_ASE_types_percentage, SpliceSeq_ASE_types_percentage)

plot_gene_names <- function(gene_counts){
  top_20_genes <- head(sort(gene_counts, decreasing = TRUE), 20)
  print(top_20_genes)
  ggplot(as.data.frame(top_20_genes), aes(x = V1, y = Freq, fill = V1)) +
    geom_col() +
    labs(title = "Histogram of Gene Frequencies ASE SpliceSeq", x = "Gene ID", y = "Number ASE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_ASE_type <- function(combined_df){
  #ASE_types_combined <- bind_rows(SplAdder_ASE_types_percentage, SpliceSeq_ASE_types_percentage)
  ggplot(combined_df, aes(x = Var1, y = Freq, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = c("SplAdder" = "aquamarine4", "SpliceSeq" = "azure3")) +
    labs(x = " ASE Type", y = "Percentage %", title = "Percentage of ASE Types for SplAdder and SpliceSeq Database") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}



get_gene_id_ensemble <- function(gene_symbols){
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Specify attributes
  attributes <- c("ensembl_gene_id", "external_gene_name")
  #gene_symbols <- gene_symbols$V1
  
  # Retrieve gene IDs
  gene_data <- getBM(attributes = attributes, filters = "external_gene_name", values = gene_symbols, mart = ensembl)
  
  # Print the results
  return(gene_data)
}