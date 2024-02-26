#!/usr/bin/Rscript

### Analysis note:
# This has not been run yet, need to get the AML count matices 


#load custom functions & packages
source("/pl/active/CSUClinHeme/users/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")

##### prepare data set #####

######### MODIFY #########

#set output name -- recommend including data and sample size
experiment <- "bm_cd34_AML_analysis_240220_log"
outName <- "allCells"

nFeature_RNA_high <- 5500
nFeature_RNA_low <- 100
percent.mt_high <- 12.5
nCount_RNA_high <- 75000
nCount_RNA_low <- 200

########## END MODIFY #########

load10x(din = "../input/", dout = "../output/s1/", outName = experiment, testQC = FALSE, removeRBC_pal = FALSE,
        nFeature_RNA_high = nFeature_RNA_high, nFeature_RNA_low = nFeature_RNA_low, percent.mt_high = percent.mt_high, 
       nCount_RNA_high = nCount_RNA_high, nCount_RNA_low = nCount_RNA_low)

#integrate the data using all of the four Seurat v5 integration methods
seu.obj <- integrateData(din = "../output/s1/", dout = "../output/s2/", outName = experiment, runAllMethods = TRUE)

# #use clustree to identify clustering parameters that appear most appropriate
# clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = experiment, 
#             test_dims = c(50,45,40), algorithm = 3, prefix = "integrated_snn_res.")

#complete data visualization
for (x in list("integrated.cca", "integrated.harmony", "integrated.joint", "integrated.rcpa")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_", x), 
                           final.dims = 45, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
                           prefix = "RNA_snn_res.", assay = "RNA", reduction = x,
                           saveRDS = F, return_obj = T, returnFeats = T,
                           features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
                          )
}

saveRDS(seu.obj, paste0("../output/s3/", outName,"_S3.rds"))
