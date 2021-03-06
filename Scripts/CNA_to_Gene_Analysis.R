# CNA_to_Gene_Analysis.R
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
setDTthreads(16)
if (!require(stringr)) {
    install.packages("stringr")
    library(stringr)
}
if (!require(parallel)) {
    install.packages("parallel")
    library(parallel)
}

# if (!require(copynumber)) {
#     BiocInstaller::biocLite("copynumber")
#     library(copynumber)
# }
# if (!require(IdeoViz)) {
#     BiocInstaller::biocLite("IdeoViz")
#     library(IdeoViz)
# }

# Get hg19 ideogram data
# ideo_table <- getIdeo("hg19")

create_segmented_genome <- function(genome = "hg19", segment_size = 10000) {
    # Create a GenomicRanges object with segment_size segments given the genome build
    if (genome == "hg19" | genome == "Grch37") {
        if (!require(BSgenome.Hsapiens.UCSC.hg19)) {
            BiocInstaller::biocLite("BSgenome.Hsapiens.UCSC.hg19")
            library(BSgenome.Hsapiens.UCSC.hg19)
        }
        cur_genome <- BSgenome.Hsapiens.UCSC.hg19
    } else if (genome == "hg38" | genome == "Grch38") {
        if (!require(BSgenome.Hsapiens.UCSC.hg38)) {
            BiocInstaller::biocLite("BSgenome.Hsapiens.UCSC.hg38")
            library(BSgenome.Hsapiens.UCSC.hg38)
        }
        cur_genome <- BSgenome.Hsapiens.UCSC.hg38
    } else {
        stop("Invalid genome; options are hg38 and hg19")
    }
    # Create chromosome name vector
    chromosome_names <- c(paste0("chr", seq(1, 22)), "chrX")
    # Get chromosome sizes from BSGenome's cur_genome
    cur_genome <- BSgenome::getSeq(x = cur_genome, names = chromosome_names)
    chr_lengths <- seqlengths(cur_genome)
    rm(cur_genome)
    # Partition/bin the genome into 10 Kb segments
    all_chroms <- list()
    for (i in 1:23) {
        cur_seq <- seq(from = 1, to = chr_lengths[i] - segment_size, by = segment_size)
        all_chroms[[i]] <- data.table(seq(from = 1, to = length(cur_seq)),
                                      rep(x = i, times = length(cur_seq)),
                                      seq(from = 1, to = chr_lengths[i] - segment_size, by = segment_size),
                                      seq(from = segment_size + 1, to = chr_lengths[i], by = segment_size))
    }
    bin_dt <- rbindlist(all_chroms)
    
    colnames(bin_dt) <- c("SeqNum", "Chromosome", "Start", "End")
    # bin_dt[Chromosome == "chr23"]$Chromosome <- "chrX"
    bin_gr <-
        makeGRangesFromDataFrame(
            df = bin_dt,
            keep.extra.columns = F,
            seqnames.field = "Chromosome",
            start.field = "Start",
            end.field = "End"
        )
    return(bin_gr)
}


bin_gr <- create_segmented_genome(genome = "hg19", segment_size = 10000)

# # Load CNA files and gene dictionary
# all_files <- list.files("SummarizedExperiments/CNA/", full.names = T)
# # Create result directory
# dir.create("CNA_GAIN_LOSS")

gene_dict <- fread("gene_dictionary.txt")
# Discard chrY
gene_dict <- gene_dict[chromosome_name %in% c(1:22, "X")]
# Change chrX and chrY to 23 and 24
gene_dict$chromosome_name[gene_dict$chromosome_name == "X"] <- "23"
gene_dict$chromosome_name <- as.integer(gene_dict$chromosome_name)

# Convert to GRanges
genome_gr <-
    makeGRangesFromDataFrame(
        df = gene_dict[, c("start_position",
                           "end_position",
                           "chromosome_name",
                           "ensembl_gene_id")],
        keep.extra.columns = T,
        ignore.strand = T,
        seqnames.field = "chromosome_name",
        start.field = "start_position",
        end.field = "end_position")


genes_in_segment <- function(segment_values, pval_less_col, pval_grtr_col, p_thresh = TRUE,
                             pval_thresh = 0.05, copy_thresh = FALSE, thresh = 0.1, gene_dict,
                             genome_granges) {
    # Find the genes in significant segments:
    if (!copy_thresh) {
        all_less_gr <- tryCatch(expr = {
            makeGRangesFromDataFrame(df = segment_values[wilcox_p_less_adj <= pval_thresh,
                                              c("seqnames", "start", "end", "wilcox_p_less_adj"), with = F],
                                     seqnames.field = "seqnames", keep.extra.columns = T)
        }, error = function(...) return(NULL))
        
        all_grtr_gr <- tryCatch(expr = {
            makeGRangesFromDataFrame(df = segment_values[wilcox_p_grtr_adj <= pval_thresh,
                                              c("seqnames", "start", "end", "wilcox_p_grtr_adj")],
                                     seqnames.field = "seqnames", keep.extra.columns = T)
        }, error = function(...) return(NULL))
    } else if (copy_thresh == T) {
        sample_cols <- colnames(all)[grepl(pattern = "GSM", x = colnames(all))]
        # Subset based on copy number thresholds (thresh difference from ploidy), count the samples
        # that fall into gain, loss and neutral category
        segment_values[, count_gain := sum(.SD >= 2 + thresh), .SDcols = sample_cols, by = "index"]
        segment_values[, count_loss := sum(.SD <= 2 - thresh), .SDcols = sample_cols, by = "index"]
        segment_values[, count_neut := length(sample_cols) - sum(count_gain, count_loss), by = "index"]
        
        all_gr <- tryCatch(expr = {
            makeGRangesFromDataFrame(df = all[, c("seqnames", "start", "end", pval_less_col, pval_grtr_col), with = F],
                                     seqnames.field = "seqnames", keep.extra.columns = T)
        }, error = function(...) return(NULL))
        
        
        all_ol <- findOverlaps(query = genome_gr, subject = all_gr, type = "any")
        # Find genome hits, classify gain and loss
        all_genes <- cbind(data.table(ensg = genome_gr[queryHits(all_ol)]$ensembl_gene_id), 
                                as.data.table(mcols(all_gr[subjectHits(all_ol), c("wilcox_p_less_adj", "wilcox_p_grtr_adj")])))
        all_genes <- unique(all_genes)
        
    }
    
    
    # Find overlaps with gene dictionary
    if (!copy_thresh) {
        loss <- NULL
        gain <- NULL
        if (!is.null(all_less_gr)) {
            all_less_ol <- findOverlaps(query = all_less_gr, subject = genome_gr, type = "any")
            # Find genome hits, classify gain and loss
            loss <- genome_gr[unique(subjectHits(all_less_ol)), ]$ensembl_gene_id
        }
        if (!is.null(all_grtr_gr)) {
            all_grtr_ol <- findOverlaps(query = all_grtr_gr, subject = genome_gr, type = "any")
            gain <- genome_gr[unique(subjectHits(all_grtr_ol)), ]$ensembl_gene_id
            # gain_pval <- all_grtr_gr[unique(queryHits(all_grtr_ol)), ]$wilcox_p_grtr_adj
        return(list(loss = loss, gain = gain))
        }
    } else {
        # count_freq <-  
    }
}

genes_in_segment <- compiler::cmpfun(genes_in_segment)
analyze_cna <- function(cur_cnv, sample_dict, study_name, destdir = getwd(), pval_thresh = 1) {
    # Find overlaps of each sample separately for less memory usage
    overlap <- findOverlaps(query = cur_cnv, subject = bin_gr,
                            type = "any", select = "all")
    final <- as.data.table(cur_cnv[queryHits(overlap)]$Segment_Mean)

    final2 <- as.data.table(bin_gr[subjectHits(overlap), ])[, c(1, 2, 3)]
    final <- cbind(final2, final, cur_cnv[queryHits(overlap)]$Sample)
    colnames(final) <- c("seqnames", "start", "end", "Segment_Mean", "sample")

    # Change shape to have samples as columns
    all <- dcast.data.table(data = final, formula = seqnames + start + end ~ sample,
                             value.var = "Segment_Mean", fun.aggregate = sum)
    # Order and hash data.table based on specified columns
    setkey(x = all, seqnames, start, end)
    
    gc()
    
    # Calculate Wilcoxon-test p-values, ensure to compare x: TP to y:NT
    wc_func <- function(cur_row) {
        less <-
            wilcox.test(x = unlist(all[cur_row, sample_dict[Types == "Tumor"]$Samples, with = F]),
                        y = unlist(all[cur_row, sample_dict[Types == "Normal"]$Samples, with = F]),
                        alternative = "less")$p.value
        greater <-
            wilcox.test(x = unlist(all[cur_row, sample_dict[Types == "Tumor"]$Samples, with = F]),
                        y = unlist(all[cur_row, sample_dict[Types == "Normal"]$Samples, with = F]),
                        alternative = "greater")$p.value
        return(data.table(less = less, greater = greater))
    }
    
    gc()
    # Compile and execute
    wc_func <- compiler::cmpfun(wc_func)
    system.time({
        wc_results <- mclapply(
            X = 1:nrow(all),
            FUN = wc_func,
            mc.preschedule = T,
            mc.cores = 8,
            mc.cleanup = T
        )
        wc_results <- rbindlist(wc_results)
        
    })
    
    # Append to all segments
    all$wilcox_p_less <- wc_results$less
    all$wilcox_p_grtr <- wc_results$greater
    # Adjust p values
    all$wilcox_p_less_adj <- p.adjust(wc_results$less, method = "fdr")
    all$wilcox_p_grtr_adj <- p.adjust(wc_results$greater, method = "fdr")
    
    all$index <- seq(1, nrow(all))
    
    # Medians and averages per tissue type
    all[, NT_median := median(unlist(.SD), na.rm = T),
                   by = "seqnames", .SDcols = sample_dict[Types == "Normal"]$Samples]
    all[, TP_median := median(unlist(.SD), na.rm = T),
                   by = "seqnames", .SDcols = sample_dict[Types == "Tumor"]$Samples]
    all[, NT_mean := mean(unlist(.SD), na.rm = T),
        by = "seqnames", .SDcols = sample_dict[Types == "Normal"]$Samples]
    all[, TP_mean := mean(unlist(.SD), na.rm = T),
        by = "seqnames", .SDcols = sample_dict[Types == "Tumor"]$Samples]
    
    # Save
    fwrite(all, paste0(destdir, "/", study_name, "_CNA_nonparam_results.txt"))
    
    # # Save the current results
    # fwrite(all, paste0(study_name, "_CNA_Genome_Overlap.csv"))
    # 
    # # Plot segment data ====
    # # merge(as.data.table(colnames(all)[4:125]), sample_dict, by.x = "V1", by.y = "Samples")
    # # ==== Plot of segment means with no thresholds ====
    # # png(filename = paste0("CNA_Plots/SegMean_Dist_No_Thresh_", proj, ".png"),
    # #     res = 220, height = 1000, width = 2000)
    # library(devtools)
    # install_github("dantaki/CNVplot")
    # library(CNVplot)
    # 
    # cur_cnv <- all[, c("seqnames", "start", "TP_mean", "NT_median", "TP_median"), with = F]
    # colnames(cur_cnv) <- c("CHROM", "POSITION", "PROBAND", "MOTHER", "FATHER")
    # CNVplot(df = cur_cnv, genome = "hg19", title = "Average Copy Number Variation")
    # 
    # plot(x = all$index, y = all$NT_mean, col = "blue",
    #                pch = 20, cex = 0.01, ylim = c(-4, 4), xaxt = "n",
    #      xlab = "", ylab = "")
    # points(x = all$index, y = all$TP_mean,
    #        col = "black", xaxt = "n", ylab = "",
    #      xlab = "", ylim = c(-4, 4), cex = 0.01)
    # for (i in 1:length(chrom_lines)) {
    #     abline(v = chrom_lines[i], col = "black", lwd = 0.1)
    # }
    # text(labels = chromosome_names,
    #      x = ((chrom_lines + chrom_lines[-1])/2)[-24],
    #      y = 0.55, srt = 45, xpd = T, adj = 0, cex = 0.6)
    # abline(h = log2(2/2), col = "red", lwd = 1)
    # 
    # title(main = paste0("Gene-wise average segment mean distribution in ", proj),
    #       ylab = "Averaged Segment Mean", mgp = c(2,0,1),
    #       xlab = paste0("Genes ordered by chromosome location. No thresholds.\nBlue: average normal (solid) | Black: average tumor (solid" ))
    # dev.off()
    
    # all <- fread("Rose_Adams_CNA_Genome_Overlap.csv")
    # Add normal and tumor segment mean medians
    # all[, normal_med := median(unlist(.SD), na.rm = T),
    #                by = c("seqnames", "start", "end"), .SDcols = all_samples[type == "NT"]$aliquot_id]
    # all[, tumor_med := median(unlist(.SD), na.rm = T),
    #     by = c("seqnames", "start", "end"), .SDcols = all_samples[type == "TP"]$aliquot_id]
    
    loss_and_gain <- genes_in_segment(segment_values = all, pval_less_col = "wilcox_p_less_adj",
                     pval_grtr_col = "wilcox_p_grtr_adj", pval_thresh = 0.01,
                     gene_dict = gene_dict, genome_granges = genome_gr)
    # Save
    saveRDS(object = loss_and_gain,
            file = paste0(destdir, "/", study_name, "_CNA_gain_loss_genes.rds"))
    
    rm(list = c("wc_results", "all_less_gr"))
    gc()
}
analyze_cna <- compiler::cmpfun(analyze_cna)

# system.time({
#     for (i in 7:length(all_files)) {
#         cat("Working on ", all_files[i], "\n")
#         analyze_cna(i)
#         gc()
#     }
# })


# ==== Find Paradoxical Genes ====
# cna_files <- list.files("CNA_GAIN_LOSS/", full.names = T)
# deseq_dea_files <- list.files("SigGenes/", pattern = "DESeq", full.names = T)
# dir.create("Paradoxical_Genes/")
# 
# idx <- 26
# find_paradoxical <- function(idx) {
#     cur <- readRDS(cna_files[idx])
#     if (is.null(cur$gain) & is.null(cur$loss)) {
#         return("No consistently amplified genes in this project")
#     }
#     # Get the project name
#     proj <- gsub(pattern = ".*(TCGA-\\w{1,4}).*",
#                  replacement = "\\1",
#                  x = cna_files[idx])
#     
#     # Find the associated file containing the differentially expressed genes for 
#     # this project
#     deseq_dea <- grepl(pattern = paste0(".*", proj, ".*"), x = deseq_dea_files)
#     
#     if (sum(deseq_dea) == 0) {
#         # The associated CNA data has no matching DEA results
#         return(paste0(proj, " has no matching DEA files"))
#     } else {
#         deseq_dea <- readRDS(deseq_dea_files[deseq_dea])
#     }
#     
#     # Find paradoxical genes
#     under_exp_parad <- cur$gain[cur$gain %in% deseq_dea$under_exp]
#     over_exp_parad <- cur$loss[cur$loss %in% deseq_dea$over_exp]
#     
#     # Save
#     saveRDS(object = list(under_exp_parad = under_exp_parad,
#                           over_exp_parad = over_exp_parad),
#             file = paste0("Paradoxical_Genes/DESeq2_paradoxical_", proj, ".rds"))
# }
# 
# find_paradoxical <- compiler::cmpfun(find_paradoxical)
# 
# mc.results <-
#     mclapply(
#         X = 1:length(cna_files),
#         FUN = find_paradoxical,
#         mc.preschedule = T,
#         mc.cores = detectCores()
#     )
# 
# 
# # ==== Plot the significant genes to contrast with threshold-setting method ====
# brca_cna <- fread("CNA_Genome_Overlaps/TCGA-BRCA_CNA_Genome_Overlap.csv")
# # Find the average segment means, by tissue type
# sample_names <- colnames(brca_cna)[-(1:3)]
# sample_split <- str_split_fixed(sample_names, "-", n = 7)
# # Get tissue type
# sample_type <- as.integer(gsub(pattern = "(\\d\\d)\\w",
#                                replacement = "\\1",
#                                x = sample_split[, 4]
# ))
# 
# # Annotate samples based on type
# # TP: Solid tumor 1
# # TM: Metastatic 6
# # NB: Blood derived normal 10
# # NT: Solid Tissue normal 11
# sample_type[sample_type == 1] <- "TP"
# sample_type[sample_type == 6] <- "TM"
# sample_type[sample_type == 10] <- "NB"
# sample_type[sample_type == 11] <- "NT"
# NT_idx <- which(sample_type == "NT")
# TP_idx <- which(sample_type == "TP")
# 
# # Order by chromosome number and start position
# data.table::setorder(x = brca_cna, seqnames, start)
# # Add index
# brca_cna$index <- 1:nrow(brca_cna)
# # Find the indices of the last gene on each chromosome
# chrom_lines <- vector(length = 23, mode = "integer")
# for (i in 1:23) {
#     chrom_lines[i] <- which(brca_cna$seqnames == i)[1]
# }
# chrom_lines <- c(chrom_lines, nrow(brca_cna))
# 
# png(filename = paste0("CNA_Plots/DESeq2_SegMean_Paradoxical_Extra", proj, ".png"),
#     res = 220, height = 1000, width = 2000)
# par(pch = 20)
# # cna_genes$color[cna_genes$color == "black"] <- "white"
# # cna_genes$color[cna_genes$color == "green"] <- "black"
# brca_cna$color <- "black"
# brca_cna[wilcox_p_less_adj <= 0.01]$color <- "red"
# brca_cna[wilcox_p_grtr_adj <= 0.01]$color <- "red"
# 
# # brca_cna$color[brca_cna$ensembl_gene_id %in% deseq_under_par] <- "red"
# # Find median by type
# brca_cna[, NT_median := median(unlist(.SD), na.rm = T),
#          by = c("start", "end", "seqnames"), .SDcols = NT_idx + 3]
# brca_cna[, TP_median := median(unlist(.SD), na.rm = T),
#          by = c("start", "end", "seqnames"), .SDcols = TP_idx + 3]
# 
# plot(x = brca_cna$index, y = brca_cna$NT_median, pch = 20,
#      col = "blue", xaxt = "n", ylab = "", cex = 0.1,
#      xlab = "", ylim = c(-.5, .5))
# points(x = brca_cna$index, y = brca_cna$TP_median, pch = 20,
#        cex = 0.1, col = "black")
# points(x = brca_cna[color == "red"]$index,
#        y = brca_cna[color == "red"]$TP_median,
#        col = "red", cex = 0.1, pch = 20)
# for (i in 1:length(chrom_lines)) {
#     abline(v = chrom_lines[i], col = "black", lwd = 0.1)
# }
# text(labels = chromosome_names,
#      x = ((chrom_lines + chrom_lines[-1])/2)[-24],
#      y = 0.55, srt = 45, xpd = T, adj = 0, cex = 0.6)
# abline(h = log2(2/2), col = "red", lwd = 0.5)
# 
# title(main = paste0("Map of aberrated regions across the genome in project TCGA-BRCA"),
#       ylab = "Medians of segment mean", mgp = c(2,0,1),
#       xlab = paste0("Mann-Whitney paired test p-values less than 0.01 are highlighted in red"))
# dev.off()
# 
