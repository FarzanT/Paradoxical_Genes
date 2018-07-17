# overlap_significance.R
# Measure the significance of overlap between two sets
source("Scripts/Package_Setup.R")
if (!require(MonteCarlo)) {
    install.packages("MonteCarlo")
    library(MonteCarlo)
}

gene_dict <- fread("gene_dictionary.txt")
# Get total number of genes, excluding the Y chromosome
gene_count <- length(unique(gene_dict[chromosome_name != "Y"]$ensembl_gene_id))


# ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# all_gene_names <- getBM(filters = "chromosome_name", attributes = ("hgnc_symbol"),
#                      values = paste0(as.character(seq(1:22)), "X"), mart = ensembl_mart)


monte_overlap_sig <- function() {
    
}

fisher_overlap_sig <- function() {
    
}

# Note that the one-tailed Fisher's exact test is the same as a hypergeometric test
fisher.test(x = ,alternative = "greater")

cur_comparison <- readRDS("Head_Neck/Head_Neck_TCGA_Comparison.rds")
cur_exp <- fread("Head_Neck/Head_Neck_limma_all_ex_results.txt")
cur_cnv <- fread("Head_Neck/Baltimore_Head_Neck_CNA_nonparam_results.txt")

cur_cn_gain_ov <- cur_comparison$cn_status$gain[wilcox_p_less_adj <= 0.05 | wilcox_p_grtr_adj <= 0.05]
cur_cn_loss_ov <- cur_comparison$cn_status$loss[wilcox_p_less_adj <= 0.05 | wilcox_p_grtr_adj <= 0.05]
cur_cn_gain <- cur_cnv[wilcox_p_grtr_adj <= 0.05]


source("Scripts/CNA_to_Gene_Analysis.R")
sub_cnv <- cur_cnv[wilcox_p_grtr_adj <= 0.05 | wilcox_p_less_adj <= 0.05]
cnv_genes <- genes_in_segment(segment_values = cur_cnv, pval_grtr_col = "wilcox_p_grtr_adj",
                 pval_less_col = "wilcox_p_less_adj", gene_dict = gene_dict,
                 genome_granges = genome_gr,
                 pval_thresh = 0.05)

gain_ov <- length(unique(cur_cn_gain$ensg))
all_cnv <- length(unique(unlist(cnv_genes)))

phyper(q = gain_ov - 1, m = all_cnv, n = gene_count - all_cnv, k = )

length(unique(cur_cn_loss$ensg))
