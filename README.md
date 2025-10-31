
getwd()
setwd("path to your folder)
getwd()
# Topics covered in this script:
# 1. Quality Control (QC) 
# 2. RMA Normalization
# 3. Pre-processing and Filtering
#INSTALL AND LOAD REQUIRED PACKAGES 
if(!requireNamespace("BiocManager", quitely = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery","affy","limma","arrayQualityMetrics"))
install.packages("dplyr")
# Load Required Libraries
library(GEOquery)             # For downloading GEO datasets
library(affy)                 # For preprocessing of Affymetrix data (RMA normalization)
library(arrayQualityMetrics)  # For generating QC reports
library(limma)
library(dplyr)                # For data manipulation
 1. Download Series Matrix Files 
# Dataset: GSE123087 (Breast Cancer validation cohort)
# Organism: Homo sapiens

# Platform: GPL25864 (Agilent-039494 SurePrint G3 Human GE v2 8x60K Microarray)

# Download series matrix
gse_data <- getGEO("GSE123087", GSEMatrix = TRUE)

# Check the structure of the GSE object
names(gse_data)

# Extract the main ExpressionSet object
eset <- gse_data[[1]]
# Extract expression data matrix (probes × samples)
expression_data <- exprs(eset)

# Extract feature data (probe annotation)
feature_data <- fData(eset)

# Extract phenotype (sample metadata)
phenotype_data <- pData(eset)

# View metadata columns
head(phenotype_data)

# Check missing values in sample annotation
sum(is.na(phenotype_data$source_name_ch1))

 2. Download and Read Raw Data (TXT files) 

# The raw data is available as TXT files inside a TAR archive
# Download from GEO (recommended to do manually if network is slow)
getGEOSuppFiles("GSE123087", baseDir = "Raw_Data", makeDirectory = TRUE)

# Untar the file
untar("Raw_Data/GSE123087/GSE123087_RAW.tar", exdir = "Raw_Data/TXT_Files")

files <- list.files("Raw_Data/TXT_Files", full.names = TRUE)
raw_data <- read.maimages(files, source = "agilent", green.only = TRUE)
summary(raw_data)
raw_data
bg_corrected <- backgroundCorrect(raw_data, method = "normexp")
normalised_data <- normalizeBetweenArrays(bg_corrected, method = "quantile")
exprs_matrix <- normalised_data$E  
dim(exprs_matrix)
eset_agilent <- ExpressionSet(assayData = exprs_matrix,
  phenotype_data = AnnotatedDataFrame(phenotype_data))


arrayQualityMetrics(expressionset = eset_agilent,
                    outdir = "Results/QC_Agilent",
                    force = TRUE,
                    do.logtransform = TRUE)
(row_median <- rowMedians(exprs_matrix)
hist(row_median, breaks = 100, freq = FALSE)
threshold <-3.5
abline(v = 3.5, col = "red", lwd = 2)
filtered_data <- exprs_matrix[row_median > threshold, ]
class(phenotype_data$source_name_ch1)
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("normal breast tissue", "breast tumor"),
                 labels = c("normal", "cancer"))
class(groups)
levels(groups)
write.csv(filtered_data, file = "Processed_Data/GSE123087_filtered_expression.csv")
saveRDS(eset_agilent, file = "Processed_Data/GSE123087_eset_agilent.rds")
save.image(file = "GSE123087filter.RData")
dim(raw_data)
