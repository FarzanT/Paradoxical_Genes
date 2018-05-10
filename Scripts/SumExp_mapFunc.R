# SumExp_mapFunc.R

# A function which takes an ExpressionSet object and returns a GRanges, or
# GRangesList object which corresponds to the genomic ranges used in the
# ExpressionSet. The rownames of the returned GRanges are used to match the
# featureNames of the ExpressionSet.

if (!require(SummarizedExperiment)) {
    install.packages("SummarizedExperiment")
    library(SummarizedExperiment)
}
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}


gene_dict <<- fread("gene_dictionary.txt")

mapFunc <- function(expSet) {
    # Extract feature names
    cur_rownames <- featureNames(expSet)
    # Find feature ranges from the gene dictionary
    data.table::merge(x = cur_rownames, y = gene_dict, by.x = )
}