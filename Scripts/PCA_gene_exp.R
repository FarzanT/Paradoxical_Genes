# PCA_gene_exp.R

if (!require(SummarizedExperiment)) {
    BiocInstaller::biocLite("SummarizedExperiment")
    library(SummarizedExperiment)
}
if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
if (!require(compiler)) {
    install.packages("compiler")
    library(compiler)
}
if (!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
}
# A function to create a PCA plot from gene expression data
pcaPlot <- function(object, projectName, intgroup = NULL, ntop = 500) {
    if (class(object) != "SummarizedExperiment") {
        object <- makeSummarizedExperimentFromExpressionSet(object)
    }
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select,]))
    percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup,
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = ":"))
    } else {
        colData(object)[[intgroup]]
    }
    d <- data.table(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = group,
        intgroup.df,
        name = colnames(object)
    )
    myPlot <- ggplot(data = d, aes(x = PC1, y = PC2, color = group)) +
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] *
                                                              100), "% variance")) + 
        ylab(paste0("PC2: ", round(percentVar[2] *
                                       100), "% variance")) + coord_fixed() +
        ggtitle(label = paste0("PCA plot of ", projectName),
                subtitle = paste0("Total number of samples used: ",
                                  ncol(object),
                                  "\nNumber of genes considered for PCA: ",
                                  ntop))+
        # Modify legend
        scale_color_discrete(name = "Tissue Type") +
        # Modify plot details
        theme(title = element_text(size = 10))
    return(myPlot)
}
pcaPlot <- compiler::cmpfun(pcaPlot)

