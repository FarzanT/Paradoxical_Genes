# GEO_Search.R
# Search GEO for useful datasets

if (!require(GEOmetadb)) {
    BiocInstaller::biocLite("GEOmetadb")
    library(GEOmetadb)
}
if (!require(GEOquery)) {
    BiocInstaller::biocLite("GEOquery")
    library(GEOquery)
}

# getSQLiteFile()

con <- dbConnect(SQLite(),'GEOmetadb.sqlite')

# Find Affymetrix SNP array 6.0 breast cancer datasets with more than 100 samples
query <- paste(
    "SELECT DISTINCT gse.title, gse.gse, gse.type",
    "FROM",
    "gsm JOIN gse_gsm ON gsm.gsm=gse_gsm.gsm",
    "JOIN gse ON gse_gsm.gse=gse.gse",
    "JOIN gse_gpl ON gse_gpl.gse=gse.gse",
    "JOIN gpl ON gse_gpl.gpl=gpl.gpl",
    "WHERE gpl.gpl = 'GPL6801'",# OR 'GPL10558' OR 'GPL570' OR 'GPL13112'",
    # "AND gse.title LIKE '%ovar%'",
    "AND gpl.organism LIKE '%Homo sapiens%'",
    "AND gse.type LIKE '%variation%'",
    "GROUP BY gse.title",
    "HAVING COUNT(*) >= 50",
    sep = " "
)
# dbGetQuery(con, "select * from gse limit 5")
(rs <- dbGetQuery(con, query))

rs$title[grep(pattern = "expression", x = rs$type, ignore.case = T)]
rs[grep(pattern = "breast", x = rs$title, ignore.case = T), c("title", "gse")]
rs$type
rs[91,]
rs[31,]

dir.create("TEMP")
# temp <- getGEO(filename = "TEMP/GSE37384_family.soft.gz", getGPL = F, AnnotGPL = F, )
temp <- getGEO(GEO = "GSE66398", getGPL = F, AnnotGPL = F)
temp <- temp$GSE66398_series_matrix.txt.gz
temp$characteristics_ch1