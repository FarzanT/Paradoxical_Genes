# Rose_Adams_Prostate_Cancer_Analysis.R

source("Scripts/Package_Setup.R")

dir.create("mySummarizedExperiments")
dir.create("limma_Results")
dir.create("Rose_Adams_Prostate")

source("limma_func.R")

# ==== Prostate Cancer ====
# Data is from GEO GSE70768
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70768

library(prostateCancerCamcap)
# library(illuminaHumanv4.db)
library(data.table)
library(SummarizedExperiment)
# raw_prostate <- getGEOSuppFiles("GSE70768")
# gunzip("GSE70768/GSE70768_non_normalized_benign.txt.gz")
# gunzip("GSE70768/GSE70768_non_normalized_tumor.txt.gz")
# untar("GSE70768/GSE70768_RAW.tar")
# fread("GSE70768/GSE70768_RAW.tar")
# benign <- fread("GSE70768/GSE70768_non_normalized_benign.txt")
# tumor <- fread("GSE70768/GSE70768_non_normalized_tumor.txt")

# TODO: USE LIMMA FOR DIFFEX, RAW COUNTS ARE NOT AVAILABLE FOR THIS DATASET

# Setup biomaRt to convert Illumina HumanHT-12 V4.0 expression beadchip IDs to
# Ensembl ENSG IDs
myMart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
all_attributes <- listAttributes(mart = myMart)
all_attributes$name[1:10]
grep(pattern = "illumina", x = all_attributes$name, ignore.case = T, value = T)

# prostate <- getGEO(GEO = "GSE70768", destdir = "GEO_Files/", GSEMatrix = T)
# prostate <- prostate$GSE70768_series_matrix.txt.gz
# Load prostate data
prostate <- prostateCancerCamcap::camcap
# Only keep Tumor and Benign data
prostate <- prostate[, prostate$Sample_Group %in% c("Tumour", "Benign")]
# Convert Illumina IDs to gene names
my_ensg <-
    biomaRt::getBM(
        mart = myMart,
        attributes = c("illumina_humanht_12_v4", "ensembl_gene_id"),
        filters = "illumina_humanht_12_v4",
        values = featureNames(prostate)
    )
my_ensg <- as.data.table(my_ensg)

# pros_ensg <-
#     as.data.table(data.frame(ENSG = unlist(
#         mget(x = featureNames(prostate), envir = illuminaHumanv4)
#     )), keep.rownames = T)
# # Remove NAs 
# pros_ensg <- pros_ensg[!is.na(pros_ensg$ENSG), ]
# Merge back with featureNames
final_names <-
    merge(x = data.table(illumina_humanht_12_v4 = featureNames(prostate)),
          y = my_ensg)

# Save for faster access
# fwrite(final_names, "illum_ensg_dict.txt")
final_names <- fread("illum_ensg_dict.txt")
# Retrieve expression data
cur_exprs <- exprs(prostate)
# Convert to data.table to merge with ENSG IDs
cur_exp_dt <- as.data.table(cur_exprs, keep.rownames = T)
colnames(cur_exp_dt)[1] <- "illumina_humanht_12_v4"
# Annotate samples by their type
colnames(cur_exp_dt)
prostate$Sample_Group


final_exprs <-
    merge(x = cur_exp_dt,
          y = final_names,
          all.x = T,
          by = "illumina_humanht_12_v4")
# Hold feature names for later addition
final_ensg <- final_exprs$ensembl_gene_id
final_illum <- final_exprs$illumina_humanht_12_v4

# Subset features that have a valid ID
# final_prostate <- prostate[featureNames(prostate) %in% final_names$rn, ]
# Remove these columns (first and last)
final_exprs <- final_exprs[, c(-1, -ncol(final_exprs)), with = F]
# Convert to data.matrix
final_exp_mat <- data.matrix(final_exprs)
rownames(final_exp_mat) <- final_ensg

# Create a SummarizedExperiment
prostate_se <-
    SummarizedExperiment(
        assays = list(RoseAdams_Prostate = final_exp_mat),
        rowData = data.table(ENSG = final_ensg, Illum_ID = final_illum),
        colData = data.table(
            Sample_ID = colnames(prostate),
            Sample_Group = prostate$Sample_Group
        )
    )
# Save 
saveRDS(prostate_se, "mySummarizedExperiments/Rose_Adams_Prostate.rds")

# ==== DEA using limma ====
library(limma)
# Load SummarizedExperiment
pros_sum_ex <- readRDS("mySummarizedExperiments/Rose_Adams_Prostate.rds")
# Use the predefined limma function to find diffferentially expressed genes
# using the given thresholds and generate an MA plot of the data
limma_func(
    cur_sum_ex = pros_sum_ex,
    cur_plot_name = "Rose-Adams Prostate Cancer",
    cur_limma_file_name = "Rose_Adams_Prostate",
    cur_min_lfc = 1,
    cur_max_pval = 0.01,
    cur_group_name = "Sample_Group",
    MA_plot = T,
    destdir = "Rose_Adams_Prostate"
)



# CNV analysis ====
# Retrieve patient CNV data from GEO
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70770
library(affy)

# BiocInstaller::biocLite("beadarray")
# bead_data <- fread("")
prostate_cnv <-
    GEOquery::getGEO(
        GEO = "GSE70770",
        destdir = "GEO_Files/",
        GSEMatrix = T,
        parseCharacteristics = T
    )

# GPL6801 (Affymetrix SNP 6.0 Array) has the samples for which we have the CEL files
prostate_beadchip <- prostate_cnv$`GSE70770-GPL6801_series_matrix.txt.gz`

featureData(prostate_cnv)
cnv_assay_data <- assayData(prostate_cnv)
phenoData(prostate_cnv)
ls(envir = cnv_assay_data)
cnv_assay_data$exprs
tmp <- fread("GEO_Files/GPL16104.soft", skip = )

dir.create("Rose_Adams_CEL_Files")

# all_files <- untar("GEO_Files/GSE70770_RAW.tar", list = T)
# cel_files <- all_files[grep(pattern = "_\\d{3,6}N|T\\.CEL\\.gz",
#                             x = all_files, ignore.case = F)]
# adj_normal <- all_files[grep(pattern = "_\\d{3,6}Adjacent\\.normal\\.CEL\\.gz",
#                              x = all_files, ignore.case = F)]
# cel_files <- c(cel_files, adj_normal)
# 
# untar(tarfile = "GEO_Files/GSE70770_RAW.tar",
#       files = cel_files, verbose = T, exdir = "GEO_Files/Rose_Adams_CEL_Files")
# cnv_batch <- affy::ReadAffy(filenames = paste0("./GEO_Files/Rose_Adams_CEL_Files/",
#                                                cel_files[1]), compress = T, verbose = T)
# 
# # Save
# saveRDS(cnv_batch, "GEO_Files/Rose_Adams_CNV_batch.rds")
# 
# slotNames(cnv_batch)
# exprs(cnv_batch)
# 
# 
# cur_idat <- limma::read.idat(idatfiles = "GEO_Files/CamCamp_Genotype_Array_1/GSE71965_RAW.tar", verbose = T)

# Everything except for Stockholm samples
all_files <- untar("GEO_Files/GSE70770_RAW.tar", list = T)
grep(pattern = ".*GSM18179\\d\\d.*", x = all_files, ignore.case = T)
"GSM1817999"
cel_files <- all_files[grep(pattern = "_\\d{3,6}N|T\\.CEL\\.gz",
                            x = all_files, ignore.case = F)]
adj_normal <- all_files[grep(pattern = "_\\d{3,6}Adjacent\\.normal\\.CEL\\.gz",
                             x = all_files, ignore.case = F)]
cel_files <- c(cel_files, adj_normal)

untar(tarfile = "GEO_Files/GSE70770_RAW.tar",
      files = cel_files, verbose = T, exdir = "GEO_Files/Rose_Adams_CEL_Files")

# Generate full paths
full_paths <- paste0("GEO_Files/Rose_Adams_CEL_Files/",
                     cel_files)
# Unzip gz CEL files, be sure to remove the '.gz' extension
for (file in full_paths) {
    gunzip(
        filename = file,
        destname = paste0(
            "Rose_Adams_Prostate/Unzipped_CEL_Files/",
            gsub(
                pattern = ".gz",
                replacement = "",
                x = basename(file))),
        remove = F)
}
# Unzipped CEL files
cel_files <- list.files("Unzipped_CEL_Files/", full.names = T)
# Use the genotype function (instead of crlmm) to allow for later analysis of
# CNA data using the CNset class. This quantile normalizes copy numbers
# The 'batch' parameter requires batch types for downstream CNV analysis; obtain
# CEL batch types from their names (T for tumor, N for normal)
batch_types <- vector(mode = "character", length = length(cel_files))
batch_types[grepl(pattern = "T\\.CEL", x = cel_files)] <- "Tumor"
# Everything else is Normal
batch_types[!grepl(pattern = "T\\.CEL", x = cel_files)] <- "Normal"

cur_cel <- crlmm::genotype(filenames = cel_files, cdfName = "genomewidesnp6",
                           batch = batch_types)

# Save
saveRDS(cur_cel, "Rose_Adams_Prostate/Rose_Adams_CNV_SNPset.rds")


# Read SNPset
cur_cel <- readRDS("Rose_Adams_Prostate/Rose_Adams_CNV_SNPset.rds")
# Change sample names to GSM IDs
sampleNames(cur_cel) <-
    gsub(pattern = "_TG.+\\.CEL",
         replacement = "",
         x = sampleNames(cur_cel))

# Find and add phenotype data from GEO files
all(sampleNames(cur_cel)%in% sampleNames(prostate_beadchip)) # TRUE

# Subset phenotype data (180 -> 129 samples)
prostate_beadchip <- prostate_beadchip[, sampleNames(cur_cel)]

phenoData(prostate_beadchip)
featureData(prostate_beadchip)

all(varLabels(cur_cel) %in% varLabels(prostate_beadchip)) # FALSE
all(fvarLabels(cur_cel) %in% fvarLabels(prostate_beadchip)) # FALSE

# Must combine the phenotype and feature data to get a complete set
phenoData(cur_cel) <- combine(phenoData(cur_cel), phenoData(prostate_beadchip))
# featureData(cur_cel) <- combine(featureData(cur_cel), featureData(prostate_beadchip))

# Save 
saveRDS(cur_cel, "Rose_Adams_Prostate/Rose_Adams_CNV_SNPset.rds")

# We now we have a CNset
class(cur_cel)
show(cur_cel)

# Perform locus- and allele-specific estimation of copy number
cur_cel.updated <- crlmmCopynumber(cur_cel)

# Re-check batch types
table(batch(cur_cel))

# Order data by chromosomal position (for faster processing)
cur_cel <- chromosomePositionOrder(cur_cel)

# Extract estimates of the raw total copy number into a matrix
raw_CNV <- rawCopynumber(cur_cel, i = seq(length = nrow(cur_cel)),
                         j = seq(length = ncol(cur_cel)))
class(raw_CNV)
# Calculate estimates of uncertainty
sds <- VanillaICE::cn_sds(raw_CNV)

# # Create an oligoClass object from CNV data, note that an 'integer matrix' of
# # raw copy number data must be provided
# oligo_CNV <- new("oligoSnpSet", cur_cel)
# oligo_CNV <- new("oligoSnpSet", copyNumber = integerMatrix(raw_CNV),
#                  featureData = featureData(cur_cel),
#                  phenoData = phenoData(cur_cel),
#                  callProbability = as.matrix(snpCallProbability(cur_cel)),
#                  call = as.matrix(calls(cur_cel)))

# Use the DNAcopy package; create a CNA object for segmentation
CNA.object <- CNA(genomdat = raw_CNV, chrom = chromosome(cur_cel),
                  maploc = position(cur_cel),
                  data.type = "logratio",
                  sampleid = sampleNames(cur_cel),
                  presorted = T)
class(CNA.object)
object.size(CNA.object) / 10e6 # ~ 200 MB

# Save
saveRDS(CNA.object, file = "Rose_Adams_Prostate/Rose_Adams_CNA_Object.rds")

# Smoothen to remove single point outliers
smoothed_CNA <- smooth.CNA(CNA.object)

# Perform circular binary segmentation (CBS) to calculated probe intensities
# across genomic regions (this will take some time)
cbs.segments <- segment(smoothed_CNA, verbose = 2)

# Save
saveRDS(cbs.segments, "Rose_Adams_Prostate/Rose_Adams_CBS.rds")


cbs.segments <- readRDS("Rose_Adams_Prostate/Rose_Adams_CBS.rds")
# cbs.segments <- readRDS("Rose_Adams_CBS.rds")
class(cbs.segments)
object.size(cbs.segments)

# Extract processed segment data
segment_data <- as.data.table(cbs.segments$output)
colnames(segment_data)
segment_data$seg.mean

# Remove rows that have NA in start, end columns or chrom
segment_data <- segment_data[!is.na(loc.start) & !is.na(loc.end) & !is.na(chrom)]
colnames(segment_data) <- c("Sample", "Chromosome", "Start", "End",
                            "Num_Probes", "Segment_Mean")
# Transform the segment means to 
# Create a GRanges object for range related analysis (e.g. finding overlaps)
cur_seqinfo <- Seqinfo(
    seqnames = as.character(seq_along(1:nrow(segment_data))),
    seqlengths = segment_data$End - segment_data$Start,
    isCircular = rep(F, nrow(segment_data)),
    genome = "hg19"
)
seg_granges <-
    makeGRangesFromDataFrame(
        df = segment_data,
        keep.extra.columns = T,
        ignore.strand = T,
        seqnames.field = "Chromosome",
        start.field = "Start",
        end.field = "End", seqinfo = cur_seqinfo
    )

# Re-annotate samples based on type (batch info from oligoClass is lost)
cel_files <- list.files("Rose_Adams_Prostate/Unzipped_CEL_Files/", full.names = T)
batch_types <- vector(mode = "character", length = length(cel_files))
batch_types[grepl(pattern = "T\\.CEL", x = cel_files)] <- "Tumor"
# Everything else is Normal
batch_types[!grepl(pattern = "T\\.CEL", x = cel_files)] <- "Normal"
samples <- gsub(pattern = ".*(GSM\\d{7})_.*", replacement = "\\1", x = cel_files)
sample_dict <- data.table(Samples = samples, Types = batch_types)

all_types <- merge(as.data.table(mcols(seg_granges)), sample_dict, by.x = "Sample", by.y = "Samples")
mcols(seg_granges) <- all_types

# Use the function defined in CNA_to_Gene_Analysis.R to measure segment means in
# 10KB segments in the hg19 genome, and then overlap significantly aberrated
# segments (found using Mann-Whitney Test) with known gene locations
# (again, in hg19)
source("Scripts/CNA_to_Gene_Analysis.R")

# Analyze CNA data and find genes in significantly aberrated regions
# RDS file will be saved as name + CNA_gain_loss 
analyze_cna(cur_cnv = seg_granges, sample_dict = sample_dict, study_name = "Rose_Adams",
            destdir = "Rose_Adams_Prostate", pval_thresh = 1)


# Differential Expression Analysis ====
all_files <- untar("GEO_Files/GSE70770_RAW.tar", list = T)
grep(pattern = ".*GSM18179\\d\\d.*", x = all_files, ignore.case = T)
"GSM1817999"
# Find idat files
idat_files <- all_files[grep(pattern = "idat",
                            x = all_files, ignore.case = F)]
dir.create("Rose_Adams_IDAT_Files")
untar(tarfile = "GEO_Files/GSE70770_RAW.tar", files = idat_files, 
      exdir = "Rose_Adams_IDAT_Files")

# Untar the bgx files
untar(tarfile = "GEO_Files/GSE70770_RAW.tar",
      files = c("GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz",
                "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz",
                "GPL20641_HumanOmni2.5M-8v1-1_B.bpm.gz",
                "GPL20641_HumanOmni25M-8v1-1_B.csv.gz"),
      exdir = "Rose_Adams_IDAT_Files")

# Unzip all files in that directory with a command
all_files <- list.files("Rose_Adams_IDAT_Files/", full.names = T)
for (file in all_files) {
    gunzip(
        filename = file,
        destname = gsub(
            pattern = ".gz",
            replacement = "",
            x = file),
        overwrite = F, remove = T)
    
}
gc()

# Read and parse HTSeq data from downloaded GEO data (takes some time)
htseq_data <-
    GEOquery::getGEO(
        GEO = "GSE70768",
        destdir = "GEO_Files/",
        GSEMatrix = T,
        parseCharacteristics = F,
        getGPL = F
    )
htseq_data <- htseq_data$GSE70768_series_matrix.txt.gz
# Add tissue type
benign <- sampleNames(htseq_data)[grepl(pattern = "benign", x = htseq_data$title)]
tumor <- sampleNames(htseq_data)[!grepl(pattern = "benign", x = htseq_data$title)]
htseq_data$`tissue_type` <- "TEMP"
htseq_data$tissue_type[grepl(pattern = "benign", x = htseq_data$title)] <- "Benign"
htseq_data$tissue_type[!grepl(pattern = "benign", x = htseq_data$title)] <- "Tumor"

# source("Scripts/PCA_gene_exp.R")
# pcaPlot(object = htseq_data, projectName = "Rose Adams", intgroup = "tissue_type", ntop = 40000)

# Note that this data is already processesed and normalized:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1817707


# tumor_htseq <-
#     GEOquery::getGEO(filename = "GEO_Files/GSE70768_non_normalized_tumor.txt.gz",
#                      destdir = "GEO_Files/",
#                      GSEMatrix = T,
#                      parseCharacteristics = T,
#                      getGPL = F
#     )
# Convert to SummarizedExperiment for use with limma
htseq_sumex <- makeSummarizedExperimentFromExpressionSet(htseq_data)
source("Scripts/limma_func.R")
limma_func(
    cur_sum_ex = htseq_sumex,
    cur_plot_name = "Rose Adams Prostate Cancer Experiment",
    cur_min_lfc = 1,
    cur_max_pval = 0.05,
    cur_group_name = "tissue_type",
    cur_limma_file_name = "Rose_Adams"
)

# Read the results
dea_results <- fread("Rose_Adams_Prostate/Rose_Adams_limma_diff_ex_results.txt")

# Find differentially expressed genes
sig_diff <- dea_results[abs(logFC) > 1 & adj.P.Val <= 0.05]

# Match Illumina probes to genes
source("Scripts/probes_to_genes.R")
cur_mapping <- probes_to_genes(probes = sig_diff$rn, type = "illumina_humanht_12_v4")
# Merge
sig_diff <- merge(sig_diff, cur_mapping, by.x = "rn", by.y = "illumina_humanht_12_v3")

# Find anti-correlated genes
cur_dea <- list(gain = sig_diff[logFC >= 1]$ensembl_gene_id, loss = sig_diff[logFC <= -1]$ensembl_gene_id)
cur_cna <- readRDS("Rose_Adams_Prostate/Rose_Adams_CNA_gain_loss_genes.rds")

source("Scripts/find_paradoxical.R")
cur_parad <- find_paradoxical(cna_list = cur_cna, dea_list = cur_dea)

# Save
saveRDS(cur_parad, "Rose_Adams_Prostate/Rose_Adams_Paradoxical_Genes.rds")

# Save LFC table with probe names replaced with ENSG IDs ====
lfc_table <- fread("Rose_Adams_Prostate/Rose_Adams_limma_all_ex_results.txt")
cur_mapping <-
    getBM(
        attributes = c(
            "illumina_humanht_12_v4",
            "ensembl_gene_id",
            "external_gene_name"
        ),
        filters = "illumina_humanht_12_v4",
        values = lfc_table$rn,
        mart = myMart
    )

lfc_table$rn <- as.character(lfc_table$rn)
lfc_table <- merge(lfc_table, cur_mapping, by.x = "rn", by.y = "illumina_humanht_12_v4")
fwrite(lfc_table, "Rose_Adams_Prostate/Rose_Adams_limma_all_ex_results.txt")

# Find the status of TCGA's paradoxical genes in this dataset ======================================
# Read TCGA-PRAD's paradoxical genes
tcga_parad <- readRDS("~/Project/Paradoxical_Genes/DESeq2_paradoxical_TCGA-PRAD.rds")

# Read the Rose-Adams diff-ex and copy number results
dea_results <- fread("Rose_Adams_Prostate/Rose_Adams_limma_all_ex_results.txt")
copy_results <- fread("Rose_Adams_Prostate/Rose_Adams_CNA_nonparam_results.txt")
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


cur_cn_status <-
    find_cn_status(
        query_gain = tcga_parad$over_exp_parad,
        query_loss = tcga_parad$under_exp_parad,
        segment_values = copy_results,
        thresh = 0.1,
        gene_dict = gene_dict,
        genome_granges = genome_gr
    )
cur_exp_status <- find_exp_status(query_gain = tcga_parad$over_exp_parad,
                                  query_loss = tcga_parad$under_exp_parad,
                                  lfc_table = dea_results, probe = TRUE,
                                  probe_type = "illumina_humanht_12_v4")

# Save
saveRDS(list(cn_status = cur_cn_status, exp_status = cur_exp_status),
        "Rose_Adams_Prostate/Rose_Adams_TCGA_Comparison.rds")
