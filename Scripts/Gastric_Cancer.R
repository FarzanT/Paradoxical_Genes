# Gastric_Cancer.R

source("Scripts/Package_Setup.R")

# Untar raw files
untar(tarfile = "Gastric/g.tar", verbose = T, exdir = "Gastric/")

cel_files <- list.files(path = "Gastric/", pattern = ".CEL.gz", full.names = T)
dir.create("Gastric/Affy_CELs")
for (file in cel_files) {
    gunzip(
        filename = file,
        destname = paste0(
            "Gastric/Affy_CELs/",
            gsub(
                pattern = ".gz",
                replacement = "",
                x = basename(file)
            )
        ),
        remove = F
    )
}

# CNV analysis =====================================================================================
# Get phenotypic data from GEO
cur_data <- getGEO(GEO = "GSE29999", destdir = "Gastric/")
cur_data <- getGEO(filename = "Gastric/GSE29999-GPL6801_series_matrix.txt.gz", AnnotGPL = T, getGPL = T)
cur_pheno <- phenoData(cur_data)

# Create a sample dictionary
sample_dict <- data.table(Samples = cur_data$geo_accession, Types = gsub(pattern = "tissue: ", replacement = "",
                                                        x = cur_data$characteristics_ch1),
                          File = paste0("Gastric/Affy_CELs//", cur_data$geo_accession, ".CEL"))

# for filename in ./*; do tail -n +1 $filename > $filename; done
# The CEL files contain a blank header; remove first line of each using command line


# Genotype using crlmm
cur_cnset <- crlmm::genotype(filenames = sample_dict$File, cdfName = "genomewidesnp6",
                             batch = sample_dict$Types, genome = "hg19", sns = sample_dict$Sample)

# Append phenotype data from GEO
phenoData(cur_cnset) <- combine(phenoData(cur_cnset), cur_pheno)

# Save
saveRDS(cur_cnset, "Gastric/Gastric_raw_CNSet.rds")

# Perform locus- and allele-specific estimation of copy number
cur_cnset.updated <- crlmmCopynumber(cur_cnset)

# Save
saveRDS(cur_cnset, "Gastric/Gastric_called_CNSet.rds")

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
saveRDS(CNA.object, file = "Gastric/Gastric_CNA_Object.rds")

# Smoothen to remove single point outliers
smoothed_CNA <- smooth.CNA(CNA.object)

# Perform circular binary segmentation (CBS) to calculated probe intensities
# across genomic regions (this will take some time)
cbs.segments <- segment(smoothed_CNA, verbose = 2)

# Save
saveRDS(cbs.segments, "Gastric/Gastric_CBS.rds")

cbs.segments <- readRDS("Gastric/Gastric_CBS.rds")

class(cbs.segments)
object.size(cbs.segments)

# Extract processed segment data
segment_data <- as.data.table(cbs.segments$output)

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
rm(list = c("raw_CNV", "smoothed_CNA", "cur_data", "cur_cnv", "cur_cnset", "cbs.segments"))
gc()
# Re-annotate samples based on type (batch info from oligoClass is lost)
all_types <- merge(as.data.table(mcols(seg_granges)), sample_dict[, c("Samples", "Types")], by.x = "Sample", by.y = "Samples")
# Update metadata columns
mcols(seg_granges) <- all_types

# Use the function defined in CNA_to_Gene_Analysis.R to measure segment means in
# 10KB segments in the hg19 genome, and then overlap significantly aberrated
# segments (found using Mann-Whitney Test) with known gene locations
# (again, in hg19)
source("Scripts/CNA_to_Gene_Analysis.R")

# Analyze CNA data and find genes in significantly aberrated regions
# RDS file will be saved as name + CNA_gain_loss
analyze_cna(cur_cnv = seg_granges, sample_dict = sample_dict,
            study_name = "Singapore_Gastric", destdir = "Gastric/", pval_thresh = 1)


# Differential Expression Analysis ===========================================================
# Read phenoData from the GEO data
source("Scripts/Package_Setup.R")
cur_expr_set <- getGEO(filename = "Gastric/GSE29999-GPL6947_series_matrix.txt.gz")
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
sampleNames(cur_expr_set) <-
    gsub(pattern = "(GSM\\d+)\\.CEL\\.gz",
         replacement = "\\1",
         x = sampleNames(cur_expr_set))


# Convert to SummarizedExperiment for use with limma
cur_sumex <- makeSummarizedExperimentFromExpressionSet(cur_expr_set)
# Add sample tissue types
sample_dict <- data.table(Samples = cur_expr_set$geo_accession,
                          Types = gsub(pattern = "tissue: ", replacement = "",
                                      x = cur_expr_set$characteristics_ch1))
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

# Save LFC table with probe names replaced with ENSG IDs ====
lfc_table <- fread("Gastric/Singapore_Gastric_limma_all_ex_results.txt")
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
fwrite(lfc_table, "Gastric/Singapore_Gastric_limma_all_ex_results.txt")

# Find anti-correlated genes
cur_dea <- list(gain = sig_diff[logFC >= 1]$ensembl_gene_id,
                loss = sig_diff[logFC <= -1]$ensembl_gene_id)
cur_cna <- readRDS("Gastric/Singapore_Gastric_CNA_gain_loss_genes.rds")

source("Scripts/find_paradoxical.R")
cur_parad <- find_paradoxical(cna_list = cur_cna, dea_list = cur_dea)

# Save
saveRDS(cur_parad, "Gastric/Singapore_Gastric_Paradoxical_Genes.rds")
