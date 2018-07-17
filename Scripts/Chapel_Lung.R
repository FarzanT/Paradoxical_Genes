# Chapel_Lung.R

source("Scripts/Package_Setup.R")

dir.create("Chapel_Lung")

# Raw data will be downloaded to the above directory from the GEO website (http)

# Untar files
dir.create("Chapel_Lung/Untarred")
untar(tarfile = "Chapel_Lung/GSE36471_RAW.tar", exdir = "Chapel_Lung/Untarred/", verbose = T)

cel_files <- list.files(path = "Chapel_Lung/Untarred/", pattern = ".CEL.gz", full.names = T)
dir.create("Chapel_Lung/Affy_CELs")
for (file in cel_files) {
    GEOquery::gunzip(
        filename = file,
        destname = paste0(
            "Chapel_Lung/Affy_CELs/",
            gsub(
                pattern = ".gz",
                replacement = "",
                x = basename(file)
            )
        ),
        remove = F
    )
}

# Get phenotypic data from GEO (or download, browser is faster)
# download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE36nnn/GSE36363/matrix/GSE36363-GPL6801_series_matrix.txt.gz",
#               destfile = "Chapel_Lung/GSE36363-GPL6801_series_matrix.txt.gz")
cur_data <- getGEO(filename = "Chapel_Lung/GSE36363-GPL6801_series_matrix.txt.gz")
cur_pheno <- phenoData(cur_data)

# Create a sample dictionary
sample_dict <- data.table(
    Samples = cur_data$geo_accession,
    Types = gsub(
        pattern = "  \\d+",
        replacement = "",
        x = cur_data$title
    ),
    File = paste0(
        "Chapel_Lung/Affy_CELs/",
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
saveRDS(cur_cnset, "Chapel_Lung/Chapel_Lung_raw_CNSet.rds")

# Perform locus- and allele-specific estimation of copy number
cur_cnset.updated <- crlmmCopynumber(cur_cnset)

# Save
saveRDS(cur_cnset, "Chapel_Lung/Chapel_Lung_CNSet.rds")

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
saveRDS(CNA.object, file = "Chapel_Lung/Chapel_Lung_CNA_Object.rds")

# Smoothen to remove single point outliers
smoothed_CNA <- smooth.CNA(CNA.object)

# Perform circular binary segmentation (CBS) to calculated probe intensities
# across genomic regions (this will take some time)
cbs.segments <- segment(smoothed_CNA, verbose = 2)

# Save
saveRDS(cbs.segments, "Chapel_Lung/Chapel_Lung_CBS.rds")

cbs.segments <- readRDS("Chapel_Lung/Chapel_Lung_CBS.rds")

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
source("Scripts/CNA_to_Gene_Analysis.R")

# Normal & Tumor types should be capitalized
sample_dict$Types <- ifelse(test = sample_dict$Types == "tumor", yes = "Tumor", no = "Normal")
# Analyze CNA data and find genes in significantly aberrated regions
# RDS file will be saved as name + CNA_gain_loss
analyze_cna(cur_cnv = seg_granges, sample_dict = sample_dict, study_name = "Chapel_Hill_Lung",
            destdir = "Chapel_Lung/", pval_thresh = 1)


# Differential Expression Analysis ===========================================================
# Read phenoData from the GEO data
cur_expr_set <- getGEO(GEO = "GSE26939")
cur_expr_set <- cur_expr_set[[1]]
cur_pheno <- phenoData(cur_expr_set)
# # Unzip CHP files
# dir.create("Gastric/Ill_Chp")
# all_ill_chp <- paste0("~/Paradoxical_Genes/",
#                       list.files("Gastric/", pattern = "\\d\\.chp\\.gz", full.names = T))
# for (file in all_ill_chp) {
#     GEOquery::gunzip(filename = file, destname = paste0("~/Paradoxical_Genes/Gastric/Ill_Chp/",
#           gsub(pattern = "\\.gz", replacement = "", x = basename(file))), remove = F, overwrite = F)
# }
# all_array_chp <- list.files("Gastric/Ill_Chp", full.names = T)
# t <- readChp(filename = all_array_chp[1])
# # Read affy btach object
# cur_affy_batch <-
#     affy::ReadAffy(
#         filenames = all_array_cel,
#         compress = T,
#         phenoData = cur_pheno,
#         verbose = T
#     )
# 
# # Use the gcrma algorithm to perform background correction, normalization and summarization
# cur_expr_set <- gcrma::gcrma(object = cur_affy_batch, fast = F)
# Update sample names
# sampleNames(cur_expr_set) <-
#     gsub(pattern = "(GSM\\d+)\\.CEL\\.gz",
#          replacement = "\\1",
#          x = sampleNames(cur_expr_set))


# Convert to SummarizedExperiment for use with limma
cur_sumex <- makeSummarizedExperimentFromExpressionSet(cur_expr_set)
# Add sample tissue types
sample_dict <- data.table(Samples = cur_expr_set$geo_accession,
                          Types = gsub(pattern = "(tumor)|", replacement = "",
                                       x = cur_expr_set$title))
sample_dict$Types <- ifelse(sample_dict$Types == "Tumor", "Tumor", "Benign")

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
    cur_plot_name = "Singapore Gastric Cancer Experiment",
    cur_min_lfc = 1,
    cur_max_pval = 0.05,
    cur_group_name = "tissue_type",
    cur_limma_file_name = "Singapore_Gastric",
    destdir = "Gastric"
)

# Read the results
dea_results <- fread("Gastric/Singapore_Gastric_limma_diff_ex_results.txt")

# Find differentially expressed genes
sig_diff <- dea_results[abs(logFC) > 1 & adj.P.Val <= 0.05]

# Match Illumina probes to genes
# Convert Affy probe IDs to Gene names, use ensembl version 90
myMart = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                 "http://Aug2017.archive.ensembl.org")

listAttributes(myMart)[1:100,]
all_filters <- listFilters(myMart)
all_filters$name[grepl(pattern = "ill", x = all_filters$description, ignore.case = T)]
cur_mapping <-
    getBM(
        attributes = c(
            "illumina_humanht_12_v4",
            "ensembl_gene_id",
            "external_gene_name"
        ),
        filters = "illumina_humanht_12_v4",
        values = sig_diff$rn,
        mart = myMart
    )
cur_mapping <- as.data.table(cur_mapping[, 1:2])

sig_diff <- merge(sig_diff, cur_mapping,
                  by.x = "rn", by.y = "illumina_humanht_12_v4")

# Find anti-correlated genes
cur_dea <- list(gain = sig_diff[logFC >= 1]$ensembl_gene_id,
                loss = sig_diff[logFC <= -1]$ensembl_gene_id)
cur_cna <- readRDS("Gastric/Singapore_Gastric_CNA_gain_loss_genes.rds")

source("Scripts/find_paradoxical.R")
cur_parad <- find_paradoxical(cna_list = cur_cna, dea_list = cur_dea)

# Save
saveRDS(cur_parad, "Gastric/Singapore_Gastric_Paradoxical_Genes.rds")