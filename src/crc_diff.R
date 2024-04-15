library(ggplot2)

# Read ASE data from SplAdder database
ASE_Spladder<- read.table("/Users/aychaserradj/Desktop/ZHAW/Praktikum/data/SplAdder_PSI_COAD_READ_short_names.tsv", header = TRUE, sep = "\t")
ASE_SpliceSeq<- read.table("/Users/aychaserradj/Desktop/ZHAW/Praktikum/data/PSI_COAD_and_READ.tsv", header = TRUE, sep = "\t")

# Get column names that end with the specific substring
# for SplAdder data base
columns_with_substring_normal <- grep("Event$|Type$|Symbol$|.11A$", names(ASE_Spladder), value = TRUE)
columns_with_substring_tumor <- grep("Event$|Type$|Symbol$|.01A$", names(ASE_Spladder), value = TRUE)

# Create a subset with tumor and normal samples
normal_samples <- subset(ASE_Spladder, select = columns_with_substring_normal)
tumor_samples <- subset(ASE_Spladder, select = columns_with_substring_tumor)

# Select normal sampls of SpliceSeq data
culumns_spliceseq_normal <- grep("symbol$|splice_type$|_Norm",  names(ASE_SpliceSeq), value = TRUE)
#  Select tumor sampls of SpliceSeq data
normal_samples_spliceseq <- subset(ASE_SpliceSeq, select = culumns_spliceseq_normal)
tumor_samples_spliceseq <- ASE_SpliceSeq[, -c(which(names(ASE_SpliceSeq) %in% culumns_spliceseq_normal))]


# Calculate mean values for each row ("na.rm" -> only on the non-missing values) 
# nan rows still preserved 
mean_values_normal <- rowMeans(normal_samples[, -(1:3)], na.rm = TRUE)
mean_values_normal_list <-  unname(as.list(mean_values_normal))
mean_values_normal <- unlist(mean_values_normal_list, use.names = FALSE)


mean_values_tumor <- rowMeans(tumor_samples[, -(1:3)], na.rm = TRUE)
mean_values_tumor_list <- unname(as.list(mean_values_tumor))
mean_values_tumor <- unlist(mean_values_tumor_list, use.names = FALSE)


# source: https://stat.ethz.ch/pipermail/r-help/2008-February/154167
perform_t_test <- function(t, n) {
  obj <- try(t.test(t, n), silent = TRUE)
  if (is(obj, "try-error")) return (NA) else return(obj$p.value)
}

# vector for p_values of the t-tests
results_p_values <- vector(mode = "numeric", length = 0)
# vector for the ASE names, to identify the most significant ASE 
ASE_names <- character(0)


# Method to perform t-test
run_t_test <- function(normal_samples, tumor_samples, db){
  results_p_values <- vector(mode = "numeric", length = 0)
  ASE_names <- character(0)
  for (i in 1:nrow(normal_samples)){
    tumor <- as.numeric(tumor_samples[i,])
    normal <- as.numeric(normal_samples[i,])
    # for spliceseq only remove first two one, SpliceAdder remove first tree one
    if (db == "spliceseq"){
      tumor_t <- tumor[-c(1,2)]
      normal_t <- normal[-c(1,2)]
    } else if (db == "spladder"){
      tumor_t <- tumor[-c(1,2,3)]
      normal_t <- normal[-c(1,2,3)]
    }
    ASE_names <- c(ASE_names, paste(normal_samples[i,][1]))
    t_test_p_value <- perform_t_test(tumor_t, normal_t)
    results_p_values <- c(results_p_values, t_test_p_value)
  }
  # Create dictionary of ASE names and p_values from t-test
  ASE_and_p_values <- setNames(results_p_values, ASE_names)
  return(ASE_and_p_values)
}

# SpliceSeq 3526 significant ASE 
# SplAdder 30859 significant ASE
eval_ASE_p_values <- function(ASE_and_p_values){
  # Empty dictionary for the significant ASE 
  ASE_p_value_significant <- c()
  # Filter to get the ASE event that are significant
  for (key in names(ASE_and_p_values)){
    if (!is.na(ASE_and_p_values[key]) && 
        !is.nan(ASE_and_p_values[key]) && 
        ASE_and_p_values[key] < 0.05){
      ASE_p_value_significant[key] <- ASE_and_p_values[key]
    }
  }
  ASE_p_value_significant_unlist <- unlist(ASE_p_value_significant)
  # sortieren funtioniert nicht so ganz
  ASE_p_value_significant_sort <- sort(ASE_p_value_significant_unlist, decreasing = FALSE)
  # TODO print ASE names
  return(ASE_p_value_significant_sort)
}


# method to create boxplot
create_box_plot <- function(normal, tumor){
  data <- list(normal=normal, tumor=tumor)
  boxplot(data, main = "Boxplot of Groups", ylab= "Values")
}

# methode ist akuell nur für SplAdder -> umbauen auch für SpliceSeq
get_ASE_genenames <- function(ASE_significant, ASE_db){
  gene_names <- c()
  ASE_significant_list <- as.list(ASE_significant)
  eval_df <- as.data.frame(ASE_significant_list)
  for (name in names(eval_df)){
    if (name %in% ASE_db$SpliceEvent){
      ASE_gene_name <- ASE_db$GeneSymbol[ASE_db$SpliceEvent == name]
      gene_names <- c(gene_names, ASE_gene_name)
  }
  }
  gene_counts <- table(gene_names)
  return(gene_counts)
}


get_AS_type <- function(ASE_significant, ASE_db){
  ASE_type <- c()
  ASE_significant_list <- as.list(ASE_significant)
  eval_df <- as.data.frame(ASE_significant_list)
  for (name in names(eval_df)){
    if (name %in% ASE_db$SpliceEvent){
      ASE_type_name <- ASE_db$SpliceType[ASE_db$SpliceEvent == name]
      ASE_type <- c(ASE_type, ASE_type_name)
    }
  }
  ASE_type_counts <- table(ASE_type)
  return(ASE_type_counts)
}


plot_gene_names <- function(gene_counts){
  top_20_genes <- head(sort(gene_counts, decreasing = TRUE), 20)
  ggplot(as.data.frame(top_20_genes), aes(x = gene_names, y = Freq, fill = gene_names)) +
    geom_col() +
    labs(title = "Histogram of Gene Frequencies", x = "Gene Names", y = "Frequency") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_ASE_type <- function(type_counts){
  ggplot(as.data.frame(type_counts), aes(x = ASE_type, y = Freq, fill = ASE_type)) +
    geom_col() +
    labs(title = "Histogram of AS Frequencies", x = "AS types", y = "Frequency") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

STR <- read.table("/Users/aychaserradj/Desktop/ZHAW/Praktikum/data/tral_and_perf_panel_meta_info_updated.tsv", header = TRUE, sep = "\t")

get_ASE_STR <- function(){
  
}
