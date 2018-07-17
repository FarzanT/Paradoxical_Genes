# Find_TCGA_status.R
# Find the copy number and gene expression status of a set of given genes in CNA and DEA datasets

source("Scripts/Package_Setup.R")

find_cn_status <- function(query_gain, query_loss, segment_values, thresh = 0.1,
                        gene_dict, genome_granges) {
    sample_cols <- colnames(segment_values)[grepl(pattern = "GSM", x = colnames(segment_values))]
    # Subset based on copy number thresholds (thresh difference from ploidy), count the samples
    # that segment_values into gain, loss and neutral category
    segment_values[, count_gain := sum(.SD >= 2 + thresh), .SDcols = sample_cols, by = "index"]
    segment_values[, count_loss := sum(.SD <= 2 - thresh), .SDcols = sample_cols, by = "index"]
    segment_values[, count_neut := length(sample_cols) - sum(count_gain, count_loss), by = "index"]
    
    all_gr <- tryCatch(expr = {
        makeGRangesFromDataFrame(df = segment_values[, c("seqnames", "start", "end",
                                                         "wilcox_p_less_adj", "wilcox_p_grtr_adj"),
                                                     with = F],
                                 seqnames.field = "seqnames", keep.extra.columns = T)
    }, error = function(...) return(NULL))
    
    
    all_ol <- findOverlaps(query = genome_gr, subject = all_gr, type = "any")
    # Find genome hits, classify gain and loss
    all_genes <- cbind(data.table(ensg = genome_gr[queryHits(all_ol)]$ensembl_gene_id), 
                       as.data.table(mcols(all_gr[subjectHits(all_ol),
                                                  c("wilcox_p_less_adj", "wilcox_p_grtr_adj")])))
    all_genes <- unique(all_genes)
    
    return(list(gain = all_genes[ensg %in% query_gain],
                loss = all_genes[ensg %in% query_loss]))
}

find_exp_status <- function(query_gain, query_loss, lfc_table, probe = TRUE, probe_type = NULL,
                            genename_col = "ENSG") {
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
        # Find overlaps based on gene name column (query is always ENSG)
        gain_overlap <- lfc_table[lfc_table[[genename_col]] %in% query_gain]
        loss_overlap <- lfc_table[lfc_table[[genename_col]] %in% query_loss]
        
        return(list(gain = gain_overlap,
                    loss = loss_overlap))
        
    }
}


comparison_table <- function(tcga_genes, cur_exp_pred, exp_merge_colname,
                             cur_cnv_pred, cnv_merge_colname) {
    # Create data.table from TCGA paradoxical gain/loss data
    all_pred <-
        rbindlist(list(
            data.table(ENSEMBL = tcga_genes$under_exp_parad, Type = "tcga_parad_loss"),
            data.table(ENSEMBL = tcga_genes$over_exp_parad, Type = "tcga_parad_gain")
        ))
    if (any(class(cur_exp_pred[[exp_merge_colname]]) == "integer")) {
        cur_exp_pred[[exp_merge_colname]] <- as.character(cur_exp_pred[[exp_merge_colname]])
    }
    # Merge with current lfc table from current dataset
    comparison <- merge(cur_exp_pred, all_pred,
                        by.x = exp_merge_colname,
                        by.y = "ENSEMBL")
    # Subset
    comparison <- comparison[, c(exp_merge_colname, "logFC", "AveExpr", "adj.P.Val", "Type"), with = F]
    
    # Find the location of TCGA genes in CNV data, return their associated copy number values:
    tcga_locations <- getBM(attributes = c("ensembl_gene_id", "start_position", "end_position", "chromosome_name"),
          filters = "ensembl_gene_id", values = all_pred$ENSEMBL, mart = myMart)
    tcga_locations$chromosome_namep[tcga_locations$chromosome_name == "X"] <- 22
    tcga_granges <- makeGRangesFromDataFrame(df = tcga_locations, keep.extra.columns = T,
                                             ignore.strand = T,
                                             seqnames.field = "chromosome_name",
                                             start.field = "start_position",
                                             end.field = "end_position")
    # Create GRanges, find overlaps, get copy number values
    cnv_granges <- makeGRangesFromDataFrame(df = cur_cnv_pred[, c("seqnames", "start", "end", "wilcox_p_less_adj",
                                                   "wilcox_p_grtr_adj", "NT_median", "TP_median",
                                                   "NT_mean", "TP_mean"), with = F],
                                            keep.extra.columns = T, ignore.strand = T,
                                            start.field = "start", end.field = "end",
                                            seqnames.field = "seqnames")
    
    cur_ov <- findOverlaps(query = tcga_granges, subject = cnv_granges, type = "any")
    
    ov_dt <- cbind(as.data.table(tcga_granges[queryHits(cur_ov), ]),
                   as.data.table(cnv_granges[subjectHits(cur_ov), ]))
    ov_genes <- unique(ov_dt[, c("ensembl_gene_id", "NT_median",
                                 "TP_median", "NT_mean", "TP_mean")])
    
    # Merge with above results
    final <- merge(comparison, ov_genes, by.x = exp_merge_colname, by.y = cnv_merge_colname, all.x = T)
    
    setorder(final)
    
    return(final)
}

comparison_table <- compiler::cmpfun(comparison_table)

gene_dict <- fread("gene_dictionary.txt")
dir.create("TCGA_Comparison")
myMart = useMart(biomart = "ensembl",
                 dataset = "hsapiens_gene_ensembl")

# HNSC comparison ====
tcga_genes <- readRDS("~/Project/Paradoxical_Genes/DESeq2_paradoxical_TCGA-HNSC.rds")
cur_exp_pred <- fread("Head_Neck/Head_Neck_limma_all_ex_results.txt")
cur_cnv_pred <- fread("Head_Neck/Baltimore_Head_Neck_CNA_nonparam_results.txt")

merger <- comparison_table(tcga_genes = tcga_genes, cur_exp_pred = cur_exp_pred, 
                           exp_merge_colname = "ENSEMBL", cur_cnv_pred = cur_cnv_pred,
                           cnv_merge_colname = "ensembl_gene_id")
merger$`Agreement` <- ifelse((merger$logFC > 0 & merger$Type == "tcga_parad_gain") | (merger$logFC < 0 & merger$Type == "tcga_parad_loss"), yes = T, no = F)
sum(merger$Agreement) / nrow(merger)

fwrite(merger, "TCGA_Comparison/Head_Neck_TCGA_comparison.txt")

# Gastric comparison ====
tcga_genes <- readRDS("~/Project/Paradoxical_Genes/DESeq2_paradoxical_TCGA-STAD.rds")
cur_exp_pred <- fread("Gastric/Singapore_Gastric_limma_all_ex_results.txt")
cur_cnv_pred <- fread("Gastric/Singapore_Gastric_CNA_nonparam_results.txt")

merger <- comparison_table(tcga_genes = tcga_genes, cur_exp_pred = cur_exp_pred,
                           exp_merge_colname = "ensembl_gene_id",
                           cur_cnv_pred = cur_cnv_pred, cnv_merge_colname = "ensembl_gene_id")
merger$`Agreement` <- ifelse(merger$logFC > 0 | merger$Type == "tcga_parad_gain", yes = T, no = F)
merger$`Agreement` <- ifelse((merger$logFC > 0 & merger$Type == "tcga_parad_gain") | (merger$logFC < 0 & merger$Type == "tcga_parad_loss"), yes = T, no = F)
sum(merger$Agreement) / nrow(merger)

# Save
fwrite(merger, "TCGA_Comparison/Gastric_TCGA_comparison.txt")


# Rose-Adams Prostate comparison ====
tcga_genes <- readRDS("~/Project/Paradoxical_Genes/DESeq2_paradoxical_TCGA-PRAD.rds")
cur_exp_pred <- fread("Rose_Adams_Prostate/Rose_Adams_limma_all_ex_results.txt")
cur_cnv_pred <- fread("Rose_Adams_Prostate/Rose_Adams_CNA_nonparam_results.txt")

merger <- comparison_table(tcga_genes = tcga_genes, cur_exp_pred = cur_exp_pred,
                           exp_merge_colname = "ensembl_gene_id",
                           cur_cnv_pred = cur_cnv_pred, cnv_merge_colname = "ensembl_gene_id")
merger$`Agreement` <- ifelse((merger$logFC > 0 & merger$Type == "tcga_parad_gain") | (merger$logFC < 0 & merger$Type == "tcga_parad_loss"), yes = T, no = F)
sum(merger$Agreement) / nrow(merger)

fwrite(merger, "TCGA_Comparison/Rose_Adams_Prostate_TCGA_comparison.txt")


# Taiwan Lung comparison ====
tcga_genes <- readRDS("~/Project/Paradoxical_Genes/DESeq2_paradoxical_TCGA-LUAD.rds")
cur_exp_pred <- fread("Lung/Taiwan_Lung_limma_all_ex_results.txt")
cur_cnv_pred <- fread("Lung/Taiwan_Lung_CNA_nonparam_results.txt")

merger <- comparison_table(tcga_genes = tcga_genes, cur_exp_pred = cur_exp_pred,
                           exp_merge_colname = "ENSEMBL",
                           cur_cnv_pred = cur_cnv_pred, cnv_merge_colname = "ensembl_gene_id")
merger$`Agreement` <- ifelse((merger$logFC > 0 & merger$Type == "tcga_parad_gain") | (merger$logFC < 0 & merger$Type == "tcga_parad_loss"), yes = T, no = F)
sum(merger$Agreement) / nrow(merger)

fwrite(merger, "TCGA_Comparison/Taiwan_Lung_TCGA_comparison.txt")


# TIGER HCC comparison ====
tcga_genes <- readRDS("~/Project/Paradoxical_Genes/DESeq2_paradoxical_TCGA-LIHC.rds")
cur_exp_pred <- fread("HCC_CCA/HCC_limma_all_ex_results.txt")
cur_cnv_pred <- fread("HCC_CCA/TIGER-HCC_CNA_nonparam_results.txt")

merger <- comparison_table(tcga_genes = tcga_genes, cur_exp_pred = cur_exp_pred,
                           exp_merge_colname = "ENSEMBL",
                           cur_cnv_pred = cur_cnv_pred, cnv_merge_colname = "ensembl_gene_id")
merger$`Agreement` <- ifelse((merger$logFC > 0 & merger$Type == "tcga_parad_gain") | (merger$logFC < 0 & merger$Type == "tcga_parad_loss"), yes = T, no = F)
sum(merger$Agreement) / nrow(merger)

fwrite(merger, "TCGA_Comparison/TIGER_HCC_TCGA_comparison.txt")

