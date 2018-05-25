# Find_TCGA_status.R
# Find the copy number and gene expression status of a set of given genes in CNA and DEA datasets

source("Scripts/Package_Setup.R")

find_cn_status <- function(query_gain, query_loss, segment_values, thresh = 0.1,
                        gene_dict, genome_granges) {
    sample_cols <- colnames(all)[grepl(pattern = "GSM", x = colnames(all))]
    # Subset based on copy number thresholds (thresh difference from ploidy), count the samples
    # that fall into gain, loss and neutral category
    all[, count_gain := sum(.SD >= 2 + thresh), .SDcols = sample_cols, by = "index"]
    all[, count_loss := sum(.SD <= 2 - thresh), .SDcols = sample_cols, by = "index"]
    all[, count_neut := length(sample_cols) - sum(count_gain, count_loss), by = "index"]
    
    all_gr <- tryCatch(expr = {
        makeGRangesFromDataFrame(df = all[, c("seqnames", "start", "end", pval_less_col, pval_grtr_col), with = F],
                                 seqnames.field = "seqnames", keep.extra.columns = T)
    }, error = function(...) return(NULL))
    
    
    all_ol <- findOverlaps(query = genome_gr, subject = all_gr, type = "any")
    # Find genome hits, classify gain and loss
    all_genes <- cbind(data.table(ensg = genome_gr[queryHits(all_ol)]$ensembl_gene_id), 
                       as.data.table(mcols(all_gr[subjectHits(all_ol), c("wilcox_p_less_adj", "wilcox_p_grtr_adj")])))
    all_genes <- unique(all_genes)
    
    return(list(all_genes[ensg %in% query_gain], all_genes[ensg %in% query_loss]))
}

find_exp_status <- function(query_gain, query_loss, lfc_table, probe = TRUE, probe_type = NULL) {
    # Given the results table (e.g. from limma), find the lfc and p-values of query genes
    if (probe == T) {
        # Convert query genes to illumina probes
        myMart = useMart(
            biomart = "ENSEMBL_MART_ENSEMBL",
            host = "grch37.ensembl.org",
            path = "/biomart/martservice",
            dataset = "hsapiens_gene_ensembl"
        )
        
        gain_illm <- as.data.table(getBM(attributes = c(probe_type, "ensembl_gene_id"),
                           filters = "ensembl_gene_id", values = query_gain, mart = myMart))
        loss_illm <- as.data.table(getBM(attributes = c(probe_type, "ensembl_gene_id"),
                           filters = "ensembl_gene_id", values = query_loss, mart = myMart))
        loss_illm <- loss_illm[loss_illm[[probe_type]] != "", ]
        gain_illm <- gain_illm[gain_illm[[probe_type]] != "", ]
        
        loss_illm[[probe_type]] <- as.character(loss_illm[[probe_type]])
        gain_illm[[probe_type]] <- as.character(gain_illm[[probe_type]])
        
        gain_merge <- merge(gain_illm, lfc_table, by.x = probe_type, by.y = "rn")
        loss_merge <- merge(loss_illm, lfc_table, by.x = probe_type, by.y = "rn")
        
        return(list(gain = gain_merge, loss = loss_merge))
    } else {
        # TODO IMPLEMENT
    }
}
