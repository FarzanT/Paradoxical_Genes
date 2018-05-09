# GEO_data

if (!require(GEOquery)) {
    install.packages("GEOquery")
    library(GEOquery)
}
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}

# Create directory to hold GEO files
dir.create("GEO_Files")
# Download a dataset of cancer sample expression data
nci_60 <- GEOquery::getGEO(GEO = "GDS4296", destdir = "GEO_Files/")

# Convert to BioConductor eSet
cur_eSet <- GDS2eSet(GDS = nci_60)
cur_breast <- cur_eSet[, cur_eSet$tissue == "breast"]
cur_breast$