# find_paradoxical.R
# Find anti-correlated genes in regards to copy number and gene expression data

if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}


# TODO REMOVE
# cna_list <- readRDS("Rose_Adams_CNA_gain_loss_genes.rds")
# dea_list <- fread("Rose_Adams_limma_diff_ex_results.txt")
# cur_mapping

find_paradoxical <- function(cna_list, dea_list) {
    # Find over-represented CNA genes in under-represented dea_list and vice versa
    parad_exp_loss <- cna_list$gain[cna_list$gain %in% dea_list$loss]
    parad_exp_gain <- cna_list$loss[cna_list$loss %in% dea_list$gain]
    
    # # Calculate p-value using Fisher's exact test (tea test)
    # matrix(data = c(length()))
    # fisher.test(x = )
    return(list(parad_exp_gain = parad_exp_gain, parad_exp_loss = parad_exp_loss))
}
