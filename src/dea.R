library("SummarizedExperiment")
library("DESeq2")
library("TCGAbiolinks")
library("biomaRt")
library("ggplot2")
library("dplyr")
# get gene expression data for Colon Adenocarcinoma tumor (TCGA-COAD) from TCGA data base
# source: https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html
download_gene_expression_data <- function(){
  
query_exp <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

GDCdownload(query_exp)

tcga_coad_summarized <- GDCprepare(query_exp, summarizedExperiment = TRUE)
coad_matrix <- assay(tcga_coad_summarized, "unstranded")
new_row_names_coad <- sub("\\..*", "", rownames(coad_matrix))
rownames(coad_matrix) <- new_row_names_coad

query_exp_read <- GDCquery(
  project = "TCGA-READ",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

GDCdownload(query_exp_read)
#GDCdownload(query_exp)

tcga_read_summarized <- GDCprepare(query_exp_read, summarizedExperiment = TRUE)
# param "unstranded" to get the count data -> necessary for differential gene expression analysis 
read_matrix <- assay(tcga_read_summarized, "unstranded")
new_row_names_read <- sub("\\..*", "", rownames(read_matrix))
rownames(read_matrix) <- new_row_names_read
}

# 41 normal samples
COAD_gene_tpm_normal <-  grepl("-11A-", colnames(coad_matrix))
COAD_gene_tpm_normal <- coad_matrix[, COAD_gene_tpm_normal]
COAD_gene_tpm_normal <- as.data.frame(COAD_gene_tpm_normal)

# 465 tumor samples
COAD_gene_tpm_tumor <-  grepl("-01A-", colnames(coad_matrix))
COAD_gene_tpm_tumor <- coad_matrix[, COAD_gene_tpm_tumor]
COAD_gene_tpm_tumor <- as.data.frame(COAD_gene_tpm_tumor)

# 10 normal samples
READ_gene_tpm_normal <-  grepl("-11A-", colnames(read_matrix))
READ_gene_tpm_normal <- read_matrix[, READ_gene_tpm_normal]
READ_gene_tpm_normal <- as.data.frame(READ_gene_tpm_normal)

# 163 tumor samples
READ_gene_tpm_tumor <-  grepl("-01A-", colnames(read_matrix))
READ_gene_tpm_tumor <- read_matrix[, READ_gene_tpm_tumor]
READ_gene_tpm_tumor <- as.data.frame(READ_gene_tpm_tumor)


normal_all <- cbind(READ_gene_tpm_normal, COAD_gene_tpm_normal)
tumor_all <- cbind(READ_gene_tpm_tumor, COAD_gene_tpm_tumor)
data_all <- cbind(normal_all, tumor_all)

data_sample_metadata <- data.frame(
  sample_id = colnames(data_all),
  condition = c(rep("normal", ncol(normal_all)), rep("tumor", ncol(tumor_all)))
)

dds <- DESeqDataSetFromMatrix(countData = data_all,
                                   colData = data_sample_metadata,
                                   design = ~ condition)
dds <- DESeq(dds)

res_all <- results(dds, name ="condition_tumor_vs_normal")

#resLFC_all <- lfcShrink(dds, coef="condition_tumor_vs_normal", type="apeglm")

summary(res_all)
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

res <- results(dds)
# Prepare data for plotting
res_df <- as.data.frame(res)
res_df$logP <- -log10(res_df$padj)
res_df <- na.omit(res_df)
res_df$significance <- ifelse(res_df$padj < 0.05, "Significant", "Not Significant")

res_df$mean <- (res_df$baseMean)
res_df$log2FoldChange <- res_df$log2FoldChange

# Remove NA values (if any)
res_df <- na.omit(res_df)
de_genes <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]

downregulated_genes <- rownames(de_genes[de_genes$log2FoldChange < 0, ])
upregulated_genes <- rownames(de_genes[de_genes$log2FoldChange > 0, ])

downregulated_genes <- as.data.frame(downregulated_genes)
upregulated_genes <- as.data.frame(upregulated_genes)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_info_upregulated <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol', 'gene_biotype'),
  filters = 'ensembl_gene_id',
  values = upregulated_genes$upregulated_genes,
  mart = ensembl
)

gene_info_downregulated <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol', "gene_biotype"),
  filters = 'ensembl_gene_id',
  values = downregulated_genes$downregulated_genes,
  mart = ensembl
)

upregulated_genes_with_symbols <- merge(upregulated_genes, gene_info_upregulated, by.x = "upregulated_genes", by.y = "ensembl_gene_id", all.x = TRUE)
downregulated_genes_with_symbols <- merge(downregulated_genes, gene_info_downregulated, by.x = "downregulated_genes", by.y = "ensembl_gene_id", all.x = TRUE)


volcano_plot <- ggplot(res_df, aes(x=log2FoldChange, y=logP, color=significance)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_color_manual(values=c("red", "blue")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )

# Print the plot
print(volcano_plot)


ma_plot <- ggplot(res_df, aes(x=mean, y=log2FoldChange, color=significance)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_x_log10() +
  scale_color_manual(values=c("red", "blue")) +
  theme_minimal() +
  labs(
    title = "MA Plot",
    x = "Mean of Normalized Counts",
    y = "Log2 Fold Change"
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )

# Print the plot
print(ma_plot)
