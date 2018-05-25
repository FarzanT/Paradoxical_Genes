# limma_func.R

if (!require(limma)) {
    BiocInstaller::biocLite("limma")
    library(limma)
}

limma_func <-
    function(cur_sum_ex,
             cur_plot_name,
             cur_min_lfc = 1,
             cur_max_pval = 0.05,
             cur_group_name,
             cur_limma_file_name,
             MA_plot = T,
             destdir = getwd()
    ) {
        
        # Convert the Sample Groups to factors (limma requirement)
        cur_sum_ex[[cur_group_name]] <- factor(cur_sum_ex[[cur_group_name]])
        
        # Create a design matrix
        cur_design <- stats::model.matrix(~ cur_sum_ex[[cur_group_name]])
        # Retrieve the first assay from SummarizedExperiment
        cur_assay <- assay(cur_sum_ex, 1)
        # Extract the counts from the SummarizedExperiment
        # count_data <- assay(cur_sum_ex, 1)
        
        # Remove genes with little or no expression values
        # count_data <- count_data[rowSums(count_data) > 1, ]
        
        # Fit the transformed object
        my_LM_fit <- lmFit(object = cur_assay, design = cur_design)
        # colnames(coef(myLMFit))
        
        # Use empirical Bayes to calculate a moderated t-statistic...
        # cat("\nRunning empirical Bayes for project", proj)
        eBayes_fit <- eBayes(fit = my_LM_fit)
        # The above object is class MArrayLM, which can be directly passed to plotMD
        
        # Select which genes are considered differentially expressed given the 
        # thresholds
        de_genes <- decideTests(eBayes_fit, method = "separate", lfc = cur_min_lfc,
                                adjust.method = "fdr", p.value = cur_max_pval)
        
        if (MA_plot == TRUE) {
            
            # Generate an averaged MA plot
            limma::plotMD(
                object = eBayes_fit,
                status = de_genes[, 2],
                main = paste0("Averaged MA plot for ", cur_plot_name, "Analysis"),
                hl.cex = 0.5,
                bg.cex = 0.5,
                sub = paste0(
                    "# Normal",
                    ": ",
                    sum(cur_sum_ex[[cur_group_name]] == "Benign"),
                    " | ",
                    "# Tumor",
                    ": ",
                    sum(cur_sum_ex[[cur_group_name]] == "Tumor"),
                    " | max p-value: ",
                    cur_max_pval,
                    " | min absolute logFC: ",
                    cur_min_lfc
                ),
                legend = F
            )
            # Add a line separating diff-ex genes
            abline(h = c(-1, 1), col = "red")
            # Add legend
            legend(
                x = "topright",
                legend = c("Over Exp", "Under Exp"),
                pch = 20,
                pt.cex = 1,
                col = c("green", "red"),
                cex = 0.8
            )
            
        }
        # Save all gene/probe information
        all_ex <- topTable(fit = eBayes_fit, coef = 2, number = Inf, sort.by = "logFC",
                           adjust.method = "BH")
        # Save the list of differentially expressed genes:
        # Select p.value and LFC thresholds, p-value adjustment methods
        # Coef = 2 refers to ???
        top_diff_ex <-
            topTable(
                fit = eBayes_fit,
                coef = 2,
                number = Inf,
                sort.by = "logFC",
                adjust.method = "BH",
                p.value = cur_max_pval,
                lfc = cur_min_lfc
            )
        
        fwrite(x = as.data.table(all_ex, keep.rownames = T),
               file = paste0(destdir, "/", cur_limma_file_name, "_limma_all_ex_results.txt"))
        # Convert to data.table and save
        fwrite(x = as.data.table(top_diff_ex, keep.rownames = T),
               file = paste0(destdir, "/", cur_limma_file_name, "_limma_diff_ex_results.txt"))
    }

# Compile
limma_func <- compiler::cmpfun(limma_func)
