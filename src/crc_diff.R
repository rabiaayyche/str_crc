
# Read ASE data from SplAdder database
ASE_Spladder<- read.table("/Users/aychaserradj/Desktop/ZHAW/Praktikum/data/SplAdder_PSI_COAD_READ_short_names.tsv", header = TRUE, sep = "\t")

# Get column names that end with the specific substring
columns_with_substring_normal <- grep("Event$|Type$|Symbol$|.11A$", names(ASE_Spladder), value = TRUE)
columns_with_substring_tumor <- grep("Event$|Type$|Symbol$|.01A$", names(ASE_Spladder), value = TRUE)

# Create a subset with tumor and normal samples
normal_samples <- subset(ASE_Spladder, select = columns_with_substring_normal)
tumor_samples <- subset(ASE_Spladder, select = columns_with_substring_tumor)


# Calculate mean values for each row ("na.rm" -> only on the non-missing values) 
# nan rows still preserved
mean_values_normal <- rowMeans(normal_samples[, -(1:3)], na.rm = TRUE)
mean_values_normal_list <-  unname(as.list(mean_values_normal))
mean_values_normal <- unlist(mean_values_normal_list, use.names = FALSE)


mean_values_tumor <- rowMeans(tumor_samples[, -(1:3)], na.rm = TRUE)
mean_values_tumor_list <- unname(as.list(mean_values_tumor))
mean_values_tumor <- unlist(mean_values_tumor_list, use.names = FALSE)

# get ASE from normal samples without any information 
nan_idx <- which(is.nan(mean_values_normal))

# Für normale samples gibt es für einige keine Informationen zu ASE -> diese ASE für die Analayse auslassen
# calculate mean values 
tumor_samples_no_nan <- subset(tumor_samples, !(seq_along(tumor_samples[, 1]) %in% nan_idx))
mean_values_tumor <- rowMeans(tumor_samples_no_nan[, -(1:3)], na.rm = TRUE)
mean_values_tumor_list <- unname(as.list(mean_values_tumor))
mean_values_tumor <- unlist(mean_values_tumor_list, use.names = FALSE)

normal_samples_no_nan <- subset(normal_samples, !(seq_along(normal_samples[, 1]) %in% nan_idx))
mean_values_normal <- rowMeans(normal_samples_no_nan[, -(1:3)], na.rm = TRUE)
mean_values_normal_list <- unname(as.list(mean_values_normal))
mean_values_normal <- unlist(mean_values_normal_list, use.names = FALSE)


# TODO: schauen welchen t-Test 
# mittelwerte berechnen für jedes ASE (normal und tumor)
perform_t_test <- function(data, alternative = "two.sided") {
  # Perform the t-test
  result <- t.test(data, alternative = alternative)
  return(result)
}
