# CHER_LOB_Study.R

# Load required packages
source("Scripts/Package_Setup.R")
# Create directory
dir.create("CHER_LOB_Dir")

# Download and parse data from GEO
cher_lob <- getGEO(
    GEO = "GSE66399",
    destdir = "CHER_LOB_Dir/",
    GSEMatrix = T,
    AnnotGPL = T,
    getGPL = T,
    parseCharacteristics = T
)
# Extract expression data
cur_exp <- cher_lob$`GSE66399-GPL570_series_matrix.txt.gz`
class(cur_exp)

# Convert ExpressionSet to SummarizedExperiment for limma compatibility
cur_exp <- makeSummarizedExperimentFromExpressionSet(cur_exp)
cur_exp$title