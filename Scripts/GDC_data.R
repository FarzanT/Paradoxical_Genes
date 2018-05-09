# GDC_data.R

if (!require(TCGAbiolinks)) {
    install.packages("TCGAbiolinks")
    library(TCGAbiolinks)
}

TCGAbiolinks::getGDCprojects()

# FM-AD -> Breast Cancer
getDataCategorySummary(project = "FM-AD")
GDCquery(
    project = "FM-AD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "HTSeq - Counts",
    legacy = F,
    access = "open"
)
TCGAbiolinks::GDCdownload(query = )