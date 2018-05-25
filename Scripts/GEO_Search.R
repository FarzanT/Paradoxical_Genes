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
geo_tables <- dbListTables(con)
geo_tables

# Find Affymetrix SNP array 6.0 breast cancer datasets with more than 100 samples
query <- paste(
    "SELECT DISTINCT gse.title, gse.gse, gse.type",
    "FROM",
    "gsm JOIN gse_gsm ON gsm.gsm=gse_gsm.gsm",
    "JOIN gse ON gse_gsm.gse=gse.gse",
    "JOIN gse_gpl ON gse_gpl.gse=gse.gse",
    "JOIN gpl ON gse_gpl.gpl=gpl.gpl",
    "WHERE gpl.gpl = 'GPL6801' OR 'GPL10558' OR 'GPL570'",
    "AND gse.title LIKE '%breast%'",
    "AND gpl.organism LIKE '%Homo sapiens%'",
    "AND gse.type LIKE '%expresssion profiling%' AND gse.type LIKE '%Genome variation%'",
    "GROUP BY gse.title",
    "HAVING COUNT(*) >= 100",
    sep = " "
)
# dbGetQuery(con, "select * from gse limit 5")
(rs <- dbGetQuery(con, query))

rs$type
rs[51,]

