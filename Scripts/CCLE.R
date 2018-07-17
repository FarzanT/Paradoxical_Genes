# CCLE.R

# Chapel_Breast.R

# Load required packages
source("Scripts/Package_Setup.R")
setwd("~/Paradoxical_Genes/")
# Create directory
dir.create("CCLE")


# Download and parse data from GEO
# The following only has the phenotype data
ccle <- getGEO(
    GEO = "GSE36138",
    destdir = "CCLE/",
    GSEMatrix = T,
    AnnotGPL = T,
    getGPL = T,
    parseCharacteristics = T
)

# Extract phenotypic data
cur_pheno <- ccle$GSE36138_series_matrix.txt.gz
table(cur_pheno$characteristics_ch1)
cur_pheno <- cur_pheno[, cur_pheno$characteristics_ch1 == "primary site: breast"]


# Untar downloaded CEL files
# (from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE76213&format=file)
dir.create("CCLE/All_CEL_Files")
all_files <- untar(tarfile = "CCLE/GSE36138_RAW.tar", list = T)

untar(tarfile = "CCLE/GSE36138_RAW.tar", exdir = "CCLE/All_CEL_Files/",
      files = basename(as.character(cur_pheno$supplementary_file)))

# Generate full paths
full_paths <- list.files("CCLE/All_CEL_Files/", full.names = T)

# Unzip gz CEL files, be sure to remove the '.gz' extension
for (file in full_paths) {
    gunzip(
        filename = file,
        destname = paste0(
            "CCLE/All_CEL_Files/",
            gsub(
                pattern = ".gz",
                replacement = "",
                x = basename(file))),
        remove = T)
}

# CNV Analysis =====================================================================================
# Load controls
gunzip(filename = "CCLE/CCLE_SNP.Arrays_2012-10-30_4_Controls.tar.gz")
dir.create("CCLE/Controls")
untar("CCLE/CCLE_SNP.Arrays_2012-10-30_4_Controls.tar", exdir = "CCLE/Controls/")
# sample_data <- fread("CCLE/CCLE_SNP.Arrays.sif_2013-12-03.txt")
sample_data <- fread("CCLE/CCLE_sample_info_file_2012-10-18.txt")
table(sample_data$`Site Primary`)
breast_samples <- sample_data[`Site Primary` == "breast"]
table(breast_samples$`Hist Subtype1`)
# Use the genotype function (instead of crlmm) to allow for later analysis of
# CNA data using the CNset class. This quantile normalizes copy numbers
# The 'batch' parameter requires batch types for downstream CNV analysis; obtain
# CEL batch types from their names (T for tumor, N for normal)
data.table(cur_pheno$characteristics_ch1.2, cur_pheno$geo_accession)
# types <- ifelse(cur_pheno$characteristics_ch1.2 == "histology subtype1: ductal_carcinoma",
#                 yes = "Tumor", no = "Normal")
cel_files <- list.files("CCLE/All_CEL_Files/", full.names = T)
control_cels <- list.files(
    "CCLE/Controls/xchip/cle/data/SNP/raw/",
    recursive = T,
    full.names = T
)
sample_dict <- data.table(
    Samples = c(cur_pheno$geo_accession, basename(control_cels)), 
    Types = c(rep("Tumor", length(cel_files)), rep("Normal", length(control_cels))),
    Files = c(paste0("CCLE/All_CEL_Files/",
                   cur_pheno$geo_accession, ".CEL"), control_cels))

# control_cels <-
#     oligo::read.celfiles(filenames = list.files(
#         "CCLE/Controls/xchip/cle/data/SNP/raw/",
#         recursive = T,
#         full.names = T
#     ))

cur_cel <- crlmm::genotype(filenames = sample_dict$Files, cdfName = "genomewidesnp6",
                           batch = sample_dict$Types, verbose = T, genome = "hg19")


# Save
saveRDS(cur_cel, "CCLE/BRCA_CNV_SNPset.rds")

# Read SNPset
cur_cnset <- readRDS("CCLE/BRCA_CNV_SNPset.rds")

# Change sample names to GSM IDs
# sampleNames(cur_cnset) <-
#     str_extract(string = sampleNames(cur_cnset), pattern = "GSM\\d+")

# Append phenotype data from GEO
# phenoData(cur_cnset) <- combine(phenoData(cur_cnset), phenoData(cur_pheno))

# Save
# saveRDS(cur_cnset, "CCLE/BRCA_raw_CNSet.rds")

# Perform locus- and allele-specific estimation of copy number
cur_cnset.updated <- crlmmCopynumber(cur_cnset)

# Save
saveRDS(cur_cnset, "CCLE/BRCA_CNSet.rds")

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
object.size(CNA.object) / 10e6 # ~420 MB

# Save
saveRDS(CNA.object, file = "CCLE/BRCA_CNA_Object.rds")

# Smoothen to remove single point outliers
smoothed_CNA <- smooth.CNA(CNA.object)

# Perform circular binary segmentation (CBS) to calculated probe intensities
# across genomic regions (this will take some time)
cbs.segments <- segment(smoothed_CNA, verbose = 2)

# Save
saveRDS(cbs.segments, "CCLE/BRCA_CBS.rds")

gc()
cbs.segments <- readRDS("CCLE/BRCA_CBS.rds")

class(cbs.segments)
object.size(cbs.segments)

# Extract processed segment data
segment_data <- as.data.table(cbs.segments$output)

# segment_data <- as.data.table(exprs(cur_pheno))
sample_dict <- data.table(Samples = sampleNames(cur_cnset),
                          Types = batch(cur_cnset))
# Fix previous GSM improper name
# segment_data$ID <- str_extract(string = segment_data$ID, pattern = "GSM\\d+")
# cbs.segments$output <- segment_data
# saveRDS(cbs.segments, "ovc_CCA/ovc_CBS.rds")

colnames(segment_data)
segment_data$ID

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
rm(list = c("raw_CNV", "smoothed_CNA", "cur_cnset", "cbs.segments", "cur_seqinfo", "cur_pheno",
            "cur_cel", "cur_cnv", "cur_cel", "CNA.object"))
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


# Analyze CNA data and find genes in significantly aberrated regions
# RDS file will be saved as name + CNA_gain_loss
detach("package:VanillaICE", unload = T)
detach("package:oligoClasses", unload = T)
unloadNamespace("oligoClasses")
unloadNamespace("crlmm")
unloadNamespace("ff")
detach("package:crlmm", unload = T)
detach("package:ff", unload = T)
gc()
analyze_cna(cur_cnv = seg_granges, sample_dict = sample_dict,
            study_name = "CCLE_Breast",
            destdir = "CCLE/", pval_thresh = 1)
gc()


# Differential Expression Analysis ===========================================================
all_files <- untar(tarfile = "CCLE/GSE36133_RAW.tar", list = T)
dir.create("CCLE/All_EXP_Files")
untar(tarfile = "CCLE/GSE36133_RAW.tar", exdir = "CCLE/All_EXP_Files/",
      files = all_files[grepl(pattern = "CEL", x = all_files)])

source("Scripts/Package_Setup.R")
setwd("~/Paradoxical_Genes/CCLE/")
gc()
brca <- getGEO(
    GEO = "GSE36133",
    destdir = "./",
    GSEMatrix = T,
    AnnotGPL = T,
    getGPL = T,
    parseCharacteristics = T
)

# Extract phenotypic data
cur_pheno <- brca$GSE36133_series_matrix.txt.gz

cur_brca <- cur_pheno$geo_accession[grepl(pattern = "breast", x = cur_pheno$characteristics_ch1)]
# Subset pheno data
cur_pheno <- cur_pheno[, cur_pheno$geo_accession %in% cur_brca]


# Add sample tissue types
sample_dict <- data.table(
    Samples = cur_pheno$geo_accession,
    Types = str_extract(string = cur_pheno$characteristics_ch1.2,
                        pattern = "(Non-Tumor)|(Tumor)"),
    File = paste0(
        "EXP_CELs/", basename(as.character(cur_pheno$supplementary_file))
    )
)

if (!require(oligo)) {
    BiocInstaller::biocLite("oligo")
    library(oligo)
}

cur_batch <-
    oligo::read.celfiles(
        filenames = sample_dict$File,
        phenoData = phenoData(cur_pheno),
        sampleNames = sample_dict$Samples,
        verbose = T
    )

# Perform robust mulit-array average; use probeset as we can analyze expression at exon levels
Sys.setenv(R_THREADS=16)
cur_expr_set <- rma(cur_batch, target = "core")

# Save
saveRDS(cur_expr_set, "brca_Expr_Set.rds")
# Use the gcrma algorithm to perform background correction, normalization and summarization
# cur_expr_set <- gcrma::gcrma(object = cur_affy_batch, fast = F)
# Update sample names
sample_dict$Types <- ifelse(sample_dict$Types == "Tumor", "Tumor", "Benign")

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
    cur_plot_name = "TIGER brca Experiment",
    cur_min_lfc = 1,
    cur_max_pval = 0.05,
    cur_group_name = "tissue_type",
    cur_limma_file_name = "brca",
    destdir = "brca_CCA"
)

# Read the results
dea_results <- fread("BC_Ovarian/brca_limma_all_ex_results.txt")

# Find differentially expressed genes
sig_diff <- dea_results[abs(logFC) > 1 & adj.P.Val <= 0.05]
sig_diff$rn <- as.character(sig_diff$rn)

# Match HuEx probes to genes
if (!require(hta20transcriptcluster.db)) {
    BiocInstaller::biocLite("hta20transcriptcluster.db")
    library(hta20transcriptcluster.db)
}

columns(hta20transcriptcluster.db)

con = hta20transcriptcluster_dbconn()
all_probes <- as.data.table(DBI::dbGetQuery(con, "select * from probes where gene_id is not null"))
sig_probes <- all_probes[probe_id %in% sig_diff$rn]
# Find the genes associated with the significant probes
sig_genes <- as.data.table(select(
    x = hta20transcriptcluster.db,
    keys = sig_probes$probe_id,
    columns = "ENSEMBL",
    keytype = "PROBEID"
))
sig_diff$rn <- as.character(sig_diff$rn)
sig_diff <- merge(sig_diff, sig_genes, by.x = "rn", by.y = "PROBEID")

# Find anti-correlated genes
cur_dea <- list(gain = sig_diff[logFC >= 1]$ENSEMBL,
                loss = sig_diff[logFC <= -1]$ENSEMBL)

cur_cna <- readRDS("brca_CCA/TIGER-brca_CNA_gain_loss_genes.rds")

source("Scripts/find_paradoxical.R")
cur_parad <- find_paradoxical(cna_list = cur_cna, dea_list = cur_dea)

# Save
saveRDS(cur_parad, "brca_CCA/brca_Paradoxical_Genes.rds")

# Save LFC table with probe names replaced with ENSG IDs ====
all_genes <- unique(as.data.table(select(
    x = hta20transcriptcluster.db,
    keys = all_probes$probe_id,
    columns = "ENSEMBL",
    keytype = "PROBEID"
)))
lfc_table <- fread("brca_CCA/ovc_limma_all_ex_results.txt")
lfc_table$rn <- as.character(lfc_table$rn)
lfc_table <- merge(lfc_table, all_genes, by.x = "rn", by.y = "PROBEID")
fwrite(lfc_table, "ovc_CCA/ovc_limma_all_ex_results.txt")
rm("lfc_table")
