
CRC_gene_tpm <- read.table("/Users/aychaserradj/Desktop/ZHAW/Praktikum/data/CRC_gene_tpm.tsv", sep = "\t", header = TRUE)


CRC_gene_name <- c("CDHR5", "UBC", "FXYD3", 
                   "RRBP1", "RPS2", "FDPS",
                   "RPL13A", "PCBP2", "LGALS4", 
                   "RPL10", "SELENBP1","ANO9", 
                   "LGALS3", "IL32", "AUP1", 
                   "RPLP0", 
                   "AP1G2", "GNB2L1", "ABHD11")
CRC_gene_id <- c("ENSG00000099834", "ENSG00000150991", "ENSG00000089356",
                 "ENSG00000125844", "ENSG00000140988", "ENSG00000160752",
                 "ENSG00000142541","ENSG00000197111", "ENSG00000171747",
                 "ENSG00000147403", "ENSG00000143416", "ENSG00000185101",
                 "ENSG00000131981", "ENSG00000008517", "ENSG00000115307",
                 "ENSG00000089157",
                 "ENSG00000213983","ENSG00000204628", "ENSG00000106077")

CRC_gene_name_id_df <- data.frame(column1 = CRC_gene_name, column2 = CRC_gene_id)

selected_rows <- CRC_gene_tpm[grep(paste0("^", paste(CRC_gene_id, collapse = "|")), CRC_gene_tpm$gene_id), ]

CRC_gene_tpm_normal <- columns_with_substring_normal <- grep("gene_id$|.11A.", names(CRC_gene_tpm), value = TRUE)
CRC_gene_tpm_normal <- subset(selected_rows, select = CRC_gene_tpm_normal)

CRC_gene_tpm_tumor <- columns_with_substring_tumor <- grep("gene_id$|.01A.", names(CRC_gene_tpm), value = TRUE)
CRC_gene_tpm_tumor <- subset(selected_rows, select = CRC_gene_tpm_tumor)


perform_t_test <- function(t, n) {
  obj <- try(t.test(t, n), silent = TRUE)
  if (is(obj, "try-error")) return (NA) else return(obj$p.value)
}

run_t_test <- function(normal_samples, tumor_samples){
  results_p_values <- vector(mode = "numeric", length = 0)
  gene_id <- character(0)
  
  results_p_values_sig <- vector(mode = "numeric", length = 0)
  gene_id_sig <- character(0)
  
  for (i in 1:nrow(normal_samples)){
    tumor <- as.numeric(tumor_samples[i,])
    normal <- as.numeric(normal_samples[i,])
    
    gene_name <- CRC_gene_tpm_normal[i,][1]

    tumor_t <- tumor[-c(1)]
    normal_t <- normal[-c(1)]
    
    t_test_p_value <- perform_t_test(tumor, normal)
    
    gene_id <- c(gene_id, gene_name)
    results_p_values <- c(results_p_values, t_test_p_value)
  }
  adjusted_p_values <- p.adjust(results_p_values, method = "BH")
  
  for (i in seq_along(adjusted_p_values)){
    if(adjusted_p_values[i] < 0.05){
      gene_id_sig <- c(gene_id_sig, gene_id[i])
      results_p_values_sig <- c(results_p_values_sig, adjusted_p_values[i])
    }
  }
  
  gene_id_p_values_sig <- setNames(results_p_values_sig, gene_id_sig)
  return(gene_id_p_values_sig)
}

