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
if (!require(limma)) {
    BiocInstaller::biocLite("limma")
    library(limma)
}



dir.create("mySummarizedExperiments")


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




results <-
    mclapply(
        X = 1:length(SumExpFiles),
        FUN = LimmaDEA,
        mc.preschedule = T,
        mc.cores = detectCores(logical = F)
    )

# Convert ExpressionSet to RangedSummarizedExperiment
# makeSummarizedExperimentFromExpressionSet(from = prostate, )
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