# Lung_Cancer.R
source("Scripts/Package_Setup.R")

dir.create("Lung")

# Attemp download from GEO
cur_data <- getGEO(GEO = "GSE33356", destdir = "Lung/", AnnotGPL = T)
# Access phenotype data
cur_pheno <- getGEO(file = "Lung/GSE33356-GPL6801_series_matrix.txt.gz")

dir.create("Lung/All_raw_files")
untar(tarfile = "Lung/GSE33356_RAW.tar", exdir = "Lung/All_raw_files", verbose = T)

# CNV Analysis =====================================================================================
all_files <- list.files("Lung/All_raw_files/", full.names = T)
all_cel_files <- all_files[grepl(pattern = "([TN]\\.CEL\\.gz)$", x = all_files)]
dir.create("Lung/Affy_CELs")
for (file in all_cel_files) {
    gunzip(
        filename = file,
        destname = paste0(
            "Lung/Affy_CELs/",
            gsub(
                pattern = ".gz",
                replacement = "",
                x = basename(file)
            )
        ),
        remove = F
    )
}
all_cel_files <- list.files(path = "Lung/Affy_CELs/", full.names = T)


batch_types <- vector(mode = "character", length = length(all_cel_files))
batch_types[grepl(pattern = "T\\.CEL", x = all_cel_files)] <- "Tumor"
# Everything else is Normal
batch_types[!grepl(pattern = "T\\.CEL", x = all_cel_files)] <- "Normal"

# Genotype
cur_cnset <- crlmm::genotype(filenames = all_cel_files, cdfName = "genomewidesnp6",
                           batch = batch_types)
# Clean sample names
sampleNames(cur_cnset) <-
    gsub(pattern = "(GSM\\d+)_.+",
         replacement = "\\1", ignore.case = T,
         x = sampleNames(cur_cnset))
# Append phenotype data from GEO
phenoData(cur_cnset) <- combine(phenoData(cur_cnset), phenoData(cur_pheno))

# Save
saveRDS(cur_cnset, "Lung/Lung_raw_CNSet.rds")

# Perform locus- and allele-specific estimation of copy number
cur_cnset.updated <- crlmmCopynumber(cur_cnset)

# Save
saveRDS(cur_cnset, "Lung/Lung_called_CNSet.rds")

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
object.size(CNA.object) / 10e6 # ~ 200 MB

# Save
saveRDS(CNA.object, file = "Lung/Lung_CNA_Object.rds")

# Smoothen to remove single point outliers
smoothed_CNA <- smooth.CNA(CNA.object)

# Perform circular binary segmentation (CBS) to calculated probe intensities
# across genomic regions (this will take some time)
cbs.segments <- segment(smoothed_CNA, verbose = 2)

# Save
saveRDS(cbs.segments, "Lung/Lung_CBS.rds")

cbs.segments <- readRDS("Lung/Lung_CBS.rds")
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
cel_files <- list.files("Lung/Affy_CELs/", full.names = T)
batch_types <- vector(mode = "character", length = length(cel_files))
batch_types[grepl(pattern = "T\\.CEL", x = cel_files)] <- "Tumor"
# Everything else is Normal
batch_types[!grepl(pattern = "T\\.CEL", x = cel_files)] <- "Normal"
samples <- gsub(pattern = ".*(GSM\\d+)_.*", replacement = "\\1", x = cel_files)
sample_dict <- data.table(Samples = samples, Types = batch_types)

all_types <- merge(as.data.table(mcols(seg_granges)), sample_dict, by.x = "Sample", by.y = "Samples")
# Update metadata columns
mcols(seg_granges) <- all_types

# Use the function defined in CNA_to_Gene_Analysis.R to measure segment means in
# 10KB segments in the hg19 genome, and then overlap significantly aberrated
# segments (found using Mann-Whitney Test) with known gene locations
# (again, in hg19)
source("Scripts/CNA_to_Gene_Analysis.R")

# Analyze CNA data and find genes in significantly aberrated regions
# RDS file will be saved as name + CNA_gain_loss 
analyze_cna(cur_cnv = seg_granges, sample_dict = sample_dict, study_name = "Taiwan_Lung",
            destdir = "Lung", pval_thresh = 1)


# Differential Expression Analysis ===========================================================
# Read phenoData from the GEO data
cur_pheno <- getGEO(filename = "Lung/GSE33356-GPL570_series_matrix.txt.gz")
cur_pheno <- phenoData(cur_pheno)
# List Affymetrix CEL files (CEL files are raw, CHP is processed)
all_array_cel <- list.files("Lung/All_raw_files/", pattern = "\\d\\.CEL\\.gz", full.names = T)
# Read affy btach object
cur_affy_batch <-
    affy::ReadAffy(
        filenames = all_array_cel,
        compress = T,
        phenoData = cur_pheno,
        verbose = T
    )

# Use the gcrma algorithm to perform background correction, normalization and summarization
cur_expr_set <- gcrma::gcrma(object = cur_affy_batch, fast = F)
# Update sample names
sampleNames(cur_expr_set) <-
    gsub(pattern = "(GSM\\d+)\\.CEL\\.gz",
         replacement = "\\1",
         x = sampleNames(cur_expr_set))


# Convert to SummarizedExperiment for use with limma
cur_sumex <- makeSummarizedExperimentFromExpressionSet(cur_expr_set)
# Add sample tissue types
sample_dict <- data.table(Samples = cur_pheno$geo_accession,
           Type = gsub(pattern = ".+([TN]$)", replacement = "\\1", x = cur_pheno$title))
sample_dict$Type <- ifelse(sample_dict$Type == "T", "Tumor", "Benign")

# batch_types <- vector(mode = "character", length = length(all_array_cel))
# batch_types[grepl(pattern = "T\\.CEL", x = all_array_cel)] <- "Tumor"
# # Everything else is Normal
# batch_types[!grepl(pattern = "T\\.CEL", x = all_array_cel)] <- "Normal"
# samples <- gsub(pattern = ".*(GSM\\d+)_.*", replacement = "\\1", x = cel_files)
# sample_dict <- data.table(Samples = samples, Types = batch_types)
# 
cur_sumex$`tissue_type` <- merge(as.data.table(colnames(cur_sumex)),
                                 sample_dict, by.x = "V1", by.y = "Samples")$Type
source("Scripts/limma_func.R")
limma_func(
    cur_sum_ex = cur_sumex,
    cur_plot_name = "Taiwan Lung Cancer Experiment",
    cur_min_lfc = 1,
    cur_max_pval = 0.05,
    cur_group_name = "tissue_type",
    cur_limma_file_name = "Taiwan_Lung",
    destdir = "Lung"
)

# Read the results
dea_results <- fread("Lung/Taiwan_Lung_limma_diff_ex_results.txt")

# Find differentially expressed genes
sig_diff <- dea_results[abs(logFC) > 1 & adj.P.Val <= 0.05]

# Match Illumina probes to genes
# Convert Affy probe IDs to Gene names, use ensembl version 90
myMart = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                 "http://Aug2017.archive.ensembl.org")

listAttributes(myMart)[1:100,]
all_filters <- listFilters(myMart)
all_filters$name[grepl(pattern = "affy", x = all_filters$description, ignore.case = T)]
cur_mapping <-
    getBM(
        attributes = c(
            "affy_hg_u133_plus_2",
            "ensembl_gene_id",
            "external_gene_name"
        ),
        filters = "affy_hg_u133_plus_2",
        values = sig_diff$rn,
        mart = myMart
    )
cur_mapping <- as.data.table(cur_mapping[, 1:2])

sig_diff <- merge(sig_diff, cur_mapping,
                        by.x = "rn", by.y = "affy_hg_u133_plus_2")

# Find anti-correlated genes
cur_dea <- list(gain = sig_diff[logFC >= 1]$ensembl_gene_id,
                loss = sig_diff[logFC <= -1]$ensembl_gene_id)
cur_cna <- readRDS("Lung/Taiwan_Lung_CNA_gain_loss_genes.rds")

source("Scripts/find_paradoxical.R")
cur_parad <- find_paradoxical(cna_list = cur_cna, dea_list = cur_dea)

# Save
saveRDS(cur_parad, "Rose_Adams_Paradoxical_Genes.rds")
