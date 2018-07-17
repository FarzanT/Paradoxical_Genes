# Baltimore_Head_Neck.R
source("Scripts/Package_Setup.R")

dir.create("Head_Neck")

cur_cnv <- getGEO("GSE33229")

# Untar files
dir.create("Affy_CNV")
dir.create("Affy_EXP")

untar("GSE33205_RAW.tar", exdir = "Affy_EXP/")
untar("GSE33229_RAW.tar", exdir = "Affy_CNV/")


# CNV Analysis =====================================================================================
cel_files <- list.files(path = "Affy_CNV/", pattern = ".CEL.gz", full.names = T)
for (file in cel_files) {
    GEOquery::gunzip(
        filename = file,
        destname = paste0(
            "Affy_CNV/",
            gsub(
                pattern = ".gz",
                replacement = "",
                x = basename(file)
            )
        ),
        remove = T
    )
}

# Get phenotypic data from GEO (or download, browser is faster)
setwd("Head_Neck/")
cur_data <- getGEO("GSE33229")
cur_data <- cur_data$GSE33229_series_matrix.txt.gz
cur_pheno <- phenoData(cur_data)

# Create a sample dictionary
sample_dict <- data.table(
    Samples = cur_data$geo_accession,
    Types = str_extract(string = cur_data$characteristics_ch1, pattern = "(Pt)|(Ctrl)"),
    File = paste0(
        "Affy_CNV/",
        gsub(
            pattern = "\\.gz",
            replacement = "",
            x = basename(as.character(cur_data$supplementary_file))
        )
    )
)

# for filename in ./*; do tail -n +1 $filename > $filename; done
# The CEL files contain a blank header; remove first line of each using command line

rm("cur_data")
gc()
# Genotype using crlmm
cur_cnset <- crlmm::genotype(filenames = sample_dict$File, cdfName = "genomewidesnp6",
                             batch = sample_dict$Types, genome = "hg19", sns = sample_dict$Sample)

# Append phenotype data from GEO
phenoData(cur_cnset) <- combine(phenoData(cur_cnset), cur_pheno)

# Save
saveRDS(cur_cnset, "Head_Neck_raw_CNSet.rds")

# Perform locus- and allele-specific estimation of copy number
cur_cnset.updated <- crlmmCopynumber(cur_cnset)

# Save
saveRDS(cur_cnset, "Head_Neck_CNSet.rds")

# Re-check batch types
table(batch(cur_cnset))

# Order data by chromosomal position (for faster processing)
cur_cnset <- chromosomePositionOrder(cur_cnset)

# Extract estimates of the raw total copy number into a matrix
raw_CNV <- rawCopynumber(cur_cnset, i = seq(length = nrow(cur_cnset)),
                         j = seq(length = ncol(cur_cnset)))
class(raw_CNV)
# Calculate estimates of uncertainty
# sds <- VanillaICE::cn_sds(object = raw_CNV)

# # Create an oligoClass object from CNV data, note that an 'integer matrix' of
# # raw copy number data must be provided
# oligo_CNV <- new("oligoSnpSet", cur_cnset)
# oligo_CNV <- new("oligoSnpSet", copyNumber = integerMatrix(raw_CNV),
#                  featureData = featureData(cur_cnset),
#                  phenoData = phenoData(cur_cnset),
#                  callProbability = as.matrix(snpCallProbability(cur_cnset)),
#                  call = as.matrix(calls(cur_cnset)))

# Use the DNAcopy package; create a CNA object for segmentation
CNA.object <- CNA(genomdat = raw_CNV, chrom = chromosome(cur_cnset),
                  maploc = position(cur_cnset),
                  data.type = "logratio",
                  sampleid = sampleNames(cur_cnset),
                  presorted = T)
class(CNA.object)
object.size(CNA.object) / 10e6 # ~ 161 MB

# Save
saveRDS(CNA.object, file = "Head_Neck_CNA_Object.rds")

# Smoothen to remove single point outliers
smoothed_CNA <- smooth.CNA(CNA.object)

# Perform circular binary segmentation (CBS) to calculated probe intensities
# across genomic regions (this will take some time)
cbs.segments <- segment(smoothed_CNA, verbose = 2)

# Save
saveRDS(cbs.segments, "Head_Neck_CBS.rds")

gc()
cbs.segments <- readRDS("Head_Neck_CBS.rds")

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
rm(list = c("raw_CNV", "smoothed_CNA", "cur_cnset", "cbs.segments"))
gc()
# Re-annotate samples based on type (batch info from oligoClass is lost)
all_types <- merge(as.data.table(mcols(seg_granges)), sample_dict[, c("Samples", "Types")],
                   by.x = "Sample", by.y = "Samples")
# Update metadata columns
mcols(seg_granges) <- all_types

# Use the function defined in CNA_to_Gene_Analysis.R to measure segment means in
# 10KB segments in the hg19 genome, and then overlap significantly aberrated
# segments (found using Mann-Whitney Test) with known gene locations
# (again, in hg19)
setwd("~/Paradoxical_Genes/")
source("Scripts/CNA_to_Gene_Analysis.R")

# Normal & Tumor types should be capitalized
sample_dict$Types <- ifelse(test = sample_dict$Types == "Pt", yes = "Tumor", no = "Normal")
# Analyze CNA data and find genes in significantly aberrated regions
# RDS file will be saved as name + CNA_gain_loss
analyze_cna(cur_cnv = seg_granges, sample_dict = sample_dict, study_name = "Baltimore_Head_Neck",
            destdir = "Head_Neck/", pval_thresh = 1)
gc()

# Differential Expression Analysis ===========================================================
source("Scripts/Package_Setup.R")
setwd("Head_Neck/")
cur_expr_set <- getGEO("GSE33205")
cur_expr_set <- cur_expr_set[[1]]
cur_pheno <- phenoData(cur_expr_set)

# Add sample tissue types
sample_dict <- data.table(
    Samples = cur_expr_set$geo_accession,
    Types = str_extract(string = cur_expr_set$characteristics_ch1, pattern = "(Pt)|(Ctrl)"),
    File = paste0(
        "Affy_EXP/", basename(as.character(cur_expr_set$supplementary_file))
        )
    )
if (!require(oligo)) {
    BiocInstaller::biocLite("oligo")
    library(oligo)
}


# Read exonset data
cur_batch <-
# Read phenoData from the GEO data
    oligo::read.celfiles(
        filenames = sample_dict$File,
        phenoData = cur_pheno,
        sampleNames = sample_dict$Samples,
        verbose = T
    )

# Perform robust mulit-array average; use probeset as we can analyze expression at exon levels
Sys.setenv(R_THREADS=16)
cur_expr_set <- rma(cur_batch, target = "core")

# Save
saveRDS(cur_expr_set, "Head_Neck_Expr_Set.rds")
# Use the gcrma algorithm to perform background correction, normalization and summarization
# cur_expr_set <- gcrma::gcrma(object = cur_affy_batch, fast = F)
# Update sample names
sample_dict$Types <- ifelse(sample_dict$Types == "Pt", "Tumor", "Benign")

# batch_types <- vector(mode = "character", length = length(all_array_cel))
# batch_types[grepl(pattern = "T\\.CEL", x = all_array_cel)] <- "Tumor"
# # Everything else is Normal
# batch_types[!grepl(pattern = "T\\.CEL", x = all_array_cel)] <- "Normal"
# samples <- gsub(pattern = ".*(GSM\\d+)_.*", replacement = "\\1", x = cel_files)
# sample_dict <- data.table(Samples = samples, Types = batch_types)
# 
cur_expr_set$`tissue_type` <- merge(as.data.table(colnames(cur_expr_set)),
                                 sample_dict, by.x = "V1", by.y = "Samples")$Type
# Convert to SummarizedExperiment for use with limma
cur_sumex <- makeSummarizedExperimentFromExpressionSet(cur_expr_set)

setwd("~/Paradoxical_Genes/")
source("Scripts/limma_func.R")
limma_func(
    cur_sum_ex = cur_sumex,
    cur_plot_name = "Baltimore Head and Neck Cancer Experiment",
    cur_min_lfc = 1,
    cur_max_pval = 0.05,
    cur_group_name = "tissue_type",
    cur_limma_file_name = "Head_Neck",
    destdir = "Head_Neck"
)

# Read the results
dea_results <- fread("Head_Neck/Head_Neck_limma_all_ex_results.txt")

# Find differentially expressed genes
sig_diff <- dea_results[abs(logFC) > 1 & adj.P.Val <= 0.05]
sig_diff$rn <- as.character(sig_diff$rn)

# Match HuEx probes to genes
if (!require(huex10sttranscriptcluster.db)) {
    BiocInstaller::biocLite("huex10sttranscriptcluster.db")
    library(huex10sttranscriptcluster.db)
}

columns(huex10sttranscriptcluster.db)

con = huex10sttranscriptcluster_dbconn()
all_probes <- as.data.table(DBI::dbGetQuery(con, "select * from probes where gene_id is not null"))
sig_probes <- all_probes[probe_id %in% sig_diff$rn]
# Find the genes associated with the significant probes
sig_genes <- as.data.table(select(
    x = huex10sttranscriptcluster.db,
    keys = sig_probes$probe_id,
    columns = c("GENENAME", "ENSEMBL"),
    keytype = "PROBEID"
))
sig_diff$rn <- as.character(sig_diff$rn)
sig_diff <- merge(sig_diff, sig_genes, by.x = "rn", by.y = "PROBEID")

# Find anti-correlated genes
cur_dea <- list(gain = sig_diff[logFC >= 1]$ENSEMBL,
                loss = sig_diff[logFC <= -1]$ENSEMBL)

cur_cna <- readRDS("Head_Neck/Baltimore_Head_Neck_CNA_gain_loss_genes.rds")

source("Scripts/find_paradoxical.R")
cur_parad <- find_paradoxical(cna_list = cur_cna, dea_list = cur_dea)

# Save
saveRDS(cur_parad, "Head_Neck/Head_Neck_Paradoxical_Genes.rds")

# Save LFC table with probe names replaced with ENSG IDs ====
all_genes <- unique(as.data.table(select(
    x = huex10sttranscriptcluster.db,
    keys = all_probes$probe_id,
    columns = "ENSEMBL",
    keytype = "PROBEID"
)))
lfc_table <- fread("Head_Neck/Head_Neck_limma_all_ex_results.txt")
lfc_table$rn <- as.character(lfc_table$rn)
lfc_table <- merge(lfc_table, all_genes, by.x = "rn", by.y = "PROBEID")
fwrite(lfc_table, "Head_Neck/Head_Neck_limma_all_ex_results.txt")
rm("lfc_table")


# Compare to TCGA ==================================================================================
source("Scripts/Find_TCGA_status.R")

# Read TCGA-PRAD and current paradoxical genes
tcga_parad <- readRDS("~/Project/Paradoxical_Genes/DESeq2_paradoxical_TCGA-HNSC.rds")
cur_parad <- readRDS("Head_Neck/Head_Neck_Paradoxical_Genes.rds")
gain_ov <- sum(cur_parad$parad_exp_gain %in% tcga_parad$over_exp_parad)
loss_ov <- sum(cur_parad$parad_exp_loss %in% tcga_parad$under_exp_parad)
parad_gain_count <- length(cur_parad$parad_exp_gain)

gain_sig <- phyper(
    q = gain_ov - 3,
    m = length(cur_parad$parad_exp_gain),
    n = 22011 - parad_gain_count,
    k = length(tcga_parad$over_exp_parad),
    lower.tail = F
)

gain_sig

# Read the Rose-Adams diff-ex and copy number results
dea_results <- fread("Head_Neck/Head_Neck_limma_all_ex_results.txt")
copy_results <- fread("Head_Neck/Baltimore_Head_Neck_CNA_nonparam_results.txt")
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
        thresh = 0.05,
        gene_dict = gene_dict,
        genome_granges = genome_gr
    )
cur_exp_status <- find_exp_status(query_gain = tcga_parad$over_exp_parad,
                                  query_loss = tcga_parad$under_exp_parad,
                                  lfc_table = dea_results, probe = FALSE,
                                  genename_col = "ENSEMBL")

# Save
saveRDS(list(cn_status = cur_cn_status, exp_status = cur_exp_status),
        "Head_Neck/Head_Neck_TCGA_Comparison.rds")

# Measure significance of overlap ==================================================================