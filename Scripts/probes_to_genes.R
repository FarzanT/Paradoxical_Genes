# probes_to_genes.R
# Map probest to Ensembl genes using biomaRt

if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
}
if (!require(biomaRt)) {
    BiocInstaller::biocLite("biomaRt")
    library(biomaRt)
}

my_mart = useMart(biomart = "ensembl",
                 dataset = "hsapiens_gene_ensembl")
# all_atts <- listAttributes(my_mart)
# all_atts$name[grep(pattern = "HumanOmni2.5", x = all_atts$description, ignore.case = T)]
# all_atts$name[grep(pattern = "HumanOmni2.5", x = all_atts$name, ignore.case = T)]
# all_atts$name[grep(pattern = "gene", x = all_atts$description, ignore.case = T)]

probes_to_genes <- function(probes, type, map_table = NULL) {
    if (is.null(map_table)) {
        # Query the biomaRt database
        cur_mapping <-
            getBM(
                attributes = c(type, "ensembl_gene_id", "external_gene_name"),
                filters = type,
                values = probes,
                mart = my_mart
            )
    } else {
        # TODO IMPLEMENT
        return(NULL)
    }
    return(cur_mapping)
}
