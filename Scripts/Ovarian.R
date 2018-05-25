# Ovarian.R
source("Scripts/Package_Setup.R")
# Load phenotype data for CNA data
cur_data <- getGEO(filename = "Ovarian/GSE19539-GPL6801_series_matrix.txt.gz")

all_cel_files <- list.files("Ovarian/All_CEL_Files/", full.names = T)
# Load CNset generated via crlmm genotyping function
cur_cels <- readRDS("Ovarian/Ovarian_SNP_Set.rds")

cur_cels