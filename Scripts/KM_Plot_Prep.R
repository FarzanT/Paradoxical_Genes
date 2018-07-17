# KM_Plot_Prep.R

if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
if (!require(biomaRt)) {
    BiocInstaller::biocLite("biomaRt")
    library(biomaRt)
}


myMart = useMart(biomart = "ensembl",
                 dataset = "hsapiens_gene_ensembl")
# Lung ====
lung <- fread("TCGA_Comparison/Taiwan_Lung_TCGA_comparison.txt")
lung_ensg <- lung[Agreement == T]$ENSEMBL

lung_ensg <- lung[Agreement == T & (logFC > 2 | logFC < -2)]$ENSEMBL

all_atts <- listAttributes(myMart)
all_atts$name[grepl(pattern = "133",
                    x = all_atts$name,
                    ignore.case = T)]
lung_probes <-
    getBM(
        attributes = c("ensembl_gene_id", "affy_hg_u133_plus_2"),
        filters = "ensembl_gene_id",
        values = lung_ensg,
        mart = myMart
    )

lung_probes <- unique(lung_probes$affy_hg_u133_plus_2)[-1]

cat(lung_probes[121:131], sep = "\n")

# Gastric ====
gastric <- fread("TCGA_Comparison/Gastric_TCGA_comparison.txt")
gastric_ensg <- gastric[Agreement == T & (logFC > 1 | logFC < -1)]$ensembl_gene_id

gastric_probes <-
    getBM(
        attributes = c("ensembl_gene_id", "affy_hg_u133_plus_2"),
        filters = "ensembl_gene_id",
        values = gastric_ensg,
        mart = myMart
    )
gastric_probes <- unique(gastric_probes$affy_hg_u133_plus_2)
cat(gastric_probes, sep = "\n")

gastric <- fread("TCGA_Comparison/")
