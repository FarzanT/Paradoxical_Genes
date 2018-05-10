# BioConductor_data.R
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
if (!require(DESeq2)) {
    install.packages("DESeq2")
    library(DESeq2)
}
if (!require(SummarizedExperiment)) {
    install.packages("SummarizedExperiment")
    library(SummarizedExperiment)
}
if (!require(GEOquery)) {
    BiocInstaller::biocLite("GEOquery")
    library(GEOquery)
}
if (!require(biomaRt)) {
    install.packages("biomaRt")
    library(biomaRt)
}
if (!require(affy)) {
    BiocInstaller::biocLite("affy")
    library(affy)
}





BiocInstaller::biocLite("curatedBreastData")
BiocInstaller::biocLite("prostateCancerCamcap")
BiocInstaller::biocLite("illuminaHumanv4.db")

# ==== Breast Cancer ====
library(curatedBreastData)

#load up datasets that are in S4 expressionSet format.
#clinical data from master clinicalTable already linked to each sample
#in these ExpressionSets in the phenoData slot.
data(curatedBreastDataExprSetList)
# Apply the post-processing steps 
all <- processExpressionSetList(exprSetList = curatedBreastDataExprSetList, numTopVarGenes = 60000)

# Merge??? 
all[[5]]


# ==== Prostate Cancer ====
# Data is from GEO GSE70768
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70768

library(prostateCancerCamcap)
library(illuminaHumanv4.db)
raw_prostate <- getGEOSuppFiles("GSE70768")
gunzip("GSE70768/GSE70768_non_normalized_benign.txt.gz")
gunzip("GSE70768/GSE70768_non_normalized_tumor.txt.gz")
untar("GSE70768/GSE70768_RAW.tar")
fread("GSE70768/GSE70768_RAW.tar")
benign <- fread("GSE70768/GSE70768_non_normalized_benign.txt")
tumor <- fread("GSE70768/GSE70768_non_normalized_tumor.txt")

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
fwrite(final_names, "illum_ensg_dict.txt")
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

# Convert ExpressionSet to RangedSummarizedExperiment
makeSummarizedExperimentFromExpressionSet(from = prostate, )
# Create a DESeq DataSet, indicate that the counts are dependent on the
# tissue type
cur_des <- DESeq2::DESeqDataSet(se = cur_se, design = ~ Sample_Group)

# We can remove the rows that have no or nearly no information about the
# amount of gene expression (across all samples)
cur_des <- cur_des[rowSums(counts(cur_des)) > 1, ]
# Remove genes with low counts across all samples

myBPPARAM <- BiocParallel::MulticoreParam(workers = detectCores())

cur_difEx <- DESeq(object = cur_des, parallel = T, BPPARAM = myBPPARAM)
prostate$Sample_Group