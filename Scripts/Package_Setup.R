# Package_Setup.R

if (!require(stringr)) {
    install.packages("stringr")
    library(stringr)
}
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
# if (!require(DESeq2)) {
#     install.packages("DESeq2")
#     library(DESeq2)
# }
# if (!require(SummarizedExperiment)) {
#     install.packages("SummarizedExperiment")
#     library(SummarizedExperiment)
# }
if (!require(GEOquery)) {
    BiocInstaller::biocLite("GEOquery")
    library(GEOquery)
}
if (!require(biomaRt)) {
    install.packages("biomaRt")
    library(biomaRt)
}
# if (!require(affy)) {
#     BiocInstaller::biocLite("affy")
#     library(affy)
# }
if (!require(limma)) {
    BiocInstaller::biocLite("limma")
    library(limma)
}
if (!require(genomewidesnp6Crlmm)) {
    BiocInstaller::biocLite("genomewidesnp6Crlmm")
    # BiocInstaller::biocLite("genomewidesnp6_cdf")
    library(genomewidesnp6Crlmm)
}
if (!require(crlmm)) {
    BiocInstaller::biocLite("crlmm")
    library(crlmm)
}
# The ff package allows better large dataset support and RAM usage
if (!require(ff)) {
    BiocInstaller::biocLite("ff")
    library(ff)
}
# VanillaICE is used to process copy number data
if (!require(VanillaICE)) {
    BiocInstaller::biocLite("VanillaICE")
    library(VanillaICE)
}
if (!require(DNAcopy)) {
    BiocInstaller::biocLite("DNAcopy")
    library(DNAcopy)
}
if (!require(GenomicRanges)) {
    BiocInstaller::biocLite("GenomicRanges")
    library(GenomicRanges)
}
# if (!require(minfi)) {
#     BiocInstaller::biocLite("minfi")
#     library(minfi)
# }
# if (!require(affxparser)) {
#     BiocInstaller::biocLite("affxparser")
#     library(affxparser)
# }
