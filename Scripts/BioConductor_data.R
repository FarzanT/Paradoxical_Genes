# BioConductor_data.R
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}

BiocInstaller::biocLite("curatedBreastData")