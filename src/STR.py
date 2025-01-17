import pandas as pd
import matplotlib.pyplot as plt

def create_coad_psi_data():
    # read SplAdder data
    read_psi = pd.read_csv('data/SplAdder_TCGA_READ_PSI.txt', sep="\t")
    read_psi["cancer_type"] = "READ"
    coad_psi = pd.read_csv("data/SplAdder_TCGA_COAD_PSI.txt", sep="\t")
    coad_psi["cancer_type"] = "COAD"

    crc_psi = coad_psi._append(read_psi)
    #crc_psi.to_csv("data/SplAdder_PSI_COAD_and_READ.tsv", sep="\t", index=False)
    crc_psi = crc_psi[['SpliceEvent', 'GeneSymbol', 'GeneID', 'SpliceType']]

    crc_psi_spliceseq = pd.read_csv("data/PSI_COAD_and_READ.tsv", sep='\t')
    crc_psi_spliceseq = crc_psi_spliceseq[['symbol', 'as_id', 'splice_type']]
    crc_psi_spliceseq = crc_psi_spliceseq.rename(columns={"symbol":"GeneSymbol"})
    
    return crc_psi, crc_psi_spliceseq


def create_STR_bed(str_panel):
    """ 
    This method creat a .bed file for the TRAL STR panel.

    Args:
        str_panel (str): Path to TRAL STR panel.
    """
    str_full_annotatted = pd.read_csv(str_panel, header=None, sep="\t")
    str_full_annotatted = str_full_annotatted.iloc[:, [0,1,2,3,5,11]]
    str_full_annotatted = str_full_annotatted.rename(columns={0: "chr", 1: "start", 2: "end", 3: "period", 5: "msa", 11: "gene_names"})
    str_full_annotatted['strand'] = '.'

    def splice_string(s):
        # Split the string by ',' and get the first part
        return s.split(',')[0]

    str_full_annotatted['msa'] = str_full_annotatted['msa'].apply(splice_string)
    str_full_annotatted['names'] = str_full_annotatted['msa'] + ';' + str_full_annotatted['gene_names']
    column_gene_msa = str_full_annotatted.pop("names")
    str_full_annotatted.insert(3, "names", column_gene_msa)
    del str_full_annotatted['gene_names']
    del str_full_annotatted['msa']
    str_full_annotatted = str_full_annotatted.drop(str_full_annotatted.index[0])

    #str_full_annotatted.to_csv('data/str_full_bed.bed', sep="\t", header=False, index=False)

# LiftOver  STR panel hg39 -> hg19
#liftOver -bedPlus=3 data/str_full_annotated_bed.bed data/hg38ToHg19.over.chain.gz data/str_full_Hg38ToHg19.bed data/unmapped.bed

def process_liftOver(liftOver_data):
    """ 
    This method process the liftOver output.

    Args:
        liftOver_data (str): Path to STR liftOver output -> rename gene_name column and delete strand column.

    Returns:
        pandas dataframe: processed STR liftOver panel.
    """
    str_liftover = pd.read_csv(liftOver_data, sep='[\t;]', names=['chr', 'start', 'end', 'msa', 'gene_name', 'period', 'strand'])
    str_liftover = str_liftover.rename(columns={"gene_name":"GeneSymbol"})
    del str_liftover['strand']
    return str_liftover


def combine_STR_splice_data(str_liftover, crc_psi_splAdder, crc_psi_spliceseq):
    """
    This method combines the STR liftOver panel and the ASE from both databases (SplAdder & SpliceSeq) based on the column "GeneSymbol".
    This makes it possible to determine in which genes STR and ASE occur.

    Args:
        str_liftover (pandas dataframe): STR liftOver panel.
        crc_psi_splAdder (pandas dataframe): ASE from splAdder database.
        crc_psi_spliceseq (pandas dataframe): ASE from SpliceSeq database.

    Returns:
        pandas dataframe: merged dataframe of STR liftOver panel with both databases.
    """
    merged_df_splAdder = pd.merge(str_liftover, crc_psi_splAdder, on='GeneSymbol', how='inner')
    merged_df_splAdder = merged_df_splAdder.drop_duplicates()
    merged_df_splAdder = merged_df_splAdder.drop_duplicates(subset=['start', 'end'], keep=False)
    #merged_df_splAdder.to_csv("data/merged_df_splAdder.tsv", sep="\t", index=False)

    merged_df_spliceSeq = pd.merge(str_liftover, crc_psi_spliceseq, on='GeneSymbol', how='inner')
    merged_df_spliceSeq = merged_df_spliceSeq.drop_duplicates()
    merged_df_spliceSeq = merged_df_spliceSeq.drop_duplicates(subset=['start', 'end'], keep=False)
    #merged_df_spliceSeq.to_csv("data/merged_df_spliceSeq.tsv", sep="\t", index=False)
    return merged_df_splAdder, merged_df_spliceSeq


def compare_spliSeq_splAdder(merged_df_splAdder, merged_df_spliceSeq, criterion1, criterion2=None):
    """
    This method compared the combined dataframes from method "def combine_STR_splice_data(str_liftover, crc_psi_splAdder, crc_psi_spliceseq)"
    If only criterion 1 is specified, the 10 most frequently occurring ones are displayed visually for each database. 
    If both criteria are specified, the number of ASEs is displayed.

    Args:
        merged_df_splAdder (pandas dataframe): ASE (splAdder database) and STR liftOver combined dataframe.
        merged_df_spliceSeq (pandas dataframe): ASE (spliceSeq database)cand STR liftOver combined dataframe.
        criterion1 (str): Criterion to compare both databases.
        criterion2 (str, optional): Criterion to compare both databases. Defaults to None.
    """
    
    if criterion2 == None:
        gene_counts_spliceSeq = merged_df_spliceSeq.groupby(criterion1).size().sort_values(ascending=False)
        gene_counts_splAdder = merged_df_splAdder.groupby(criterion1).size().sort_values(ascending=False)
        
        # Select the 10 most frequent ones
        top_10_entries_spliceSeq = gene_counts_spliceSeq.head(10)
        top_10_entries_splAdder = gene_counts_splAdder.head(10)
    else:
        gene_counts_spliceSeq = merged_df_spliceSeq.groupby(criterion1).size().sort_values(ascending=False)
        gene_counts_splAdder = merged_df_splAdder.groupby(criterion2).size().sort_values(ascending=False)   

    plt.figure(figsize=(10, 5))

    plt.subplot(1, 2, 1)
    if criterion2 == None:
        top_10_entries_spliceSeq.plot(kind='bar', color='navy')
    else:
         gene_counts_spliceSeq.plot(kind='bar', color='navy')   
    plt.title(f'Top 10 Most Frequent {criterion1} (SpliceSeq)')
    plt.xlabel(criterion1)
    plt.ylabel('Number STR')
    plt.xticks(rotation=45)


    plt.subplot(1, 2, 2)
    if criterion2 == None:
        top_10_entries_splAdder.plot(kind='bar', color='lightsteelblue')
    else:
        gene_counts_splAdder.plot(kind='bar', color='lightsteelblue')
    plt.title(f'Top 10 Most Frequent {criterion1} (SplAdder)')
    plt.xlabel(criterion1)
    plt.ylabel('Number STR')
    plt.xticks(rotation=45)
    

    plt.tight_layout()
    #plt.show()
    plt.savefig(f"data/comparison_{criterion1}_spliceSeq_splAdder.pdf")
    
    
if __name__=="__main__": 
    crc_psi_splAdder, crc_psi_spliceSeq = create_coad_psi_data()
    create_STR_bed('data/tral_and_perf_panel_meta_info_updated.tsv')
    str_liftOver = process_liftOver("data/str_full_Hg38ToHg19.bed")
    merged_df_splAdder, merged_df_spliceSeq = combine_STR_splice_data(str_liftOver, crc_psi_splAdder, crc_psi_spliceSeq)
    #compare_spliSeq_splAdder(merged_df_splAdder, merged_df_spliceSeq, 'GeneSymbol')
    #compare_spliSeq_splAdder(merged_df_splAdder, merged_df_spliceSeq, 'splice_type', 'SpliceType')