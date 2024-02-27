#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/CSUClinHeme/users/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")

#load in processed data from https://doi.org/10.1158/2159-8290.CD-22-1297 -- reprocessed
seu.obj <- readRDS("../output/s3/ssp_mlsc_res0.6_dims45_dist0.1_neigh10_S3.rds")
outName <- "human_k9_comp"

#inspect data for proper import -- looks good
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated", 
              group.by = "scArches_Cluster",
              pt.size = 0.1,
              label = T,
              label.box = T,
              repel = F,
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_scArches_UMAP.png"), width = 7, height = 7)


#plot m-LSC score
p <- FeaturePlot(seu.obj, features = "m.LSC.score", reduction = "umap.integrated") 
p <- formatUMAP(p) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_mLSC_score_UMAP.png"), width = 7, height = 7)


#plot m-LSC score by cluster
VlnPlot(seu.obj, features = c("m.LSC.score","p.LSC.score","CD276"), group.by = "scArches_Cluster", stack = T, flip = T, fill.by = "ident")
ggsave(paste0("../output/", outName, "/", outName, "_mLSC_score_viln.png"), width = 7, height = 5)

### Key feature plots
features <- c("PTPRC","CD3E","ANPEP", 
                "DLA-DRA","CSF3R","S100A12", 
                "CD68","FLT3","FCER1A", 
                "EPCAM","COL1A1","CD34",
                "COL1A2","MS4A1","TOP2A")
p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, title.size = 14, features = features, order = F, legJust = "top", reduction = "umap.integrated") 
ggsave(paste0("../output/", outName, "/", outName, "_featPlots.png"), width = 9, height = 15)

### Key feature plots
features <- c("CD276")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 1, title.size = 14, features = features, order = T, legJust = "top", reduction = "umap.integrated") 
ggsave(paste0("../output/", outName, "/", outName, "_featPlots.png"), width = 4, height = 4)


### Integrate cross-species

#read in processed k9 data
seu.obj.k9 <- readRDS("../output/s3/can_aml_res0.6_dims45_dist0.1_neigh10_S3.rds")

#split then merge objects
message(paste0(Sys.time(), " INFO: splitting data from k9 and human."))
seu.list <- c(SplitObject(seu.obj.k9, split.by = "orig.ident"), SplitObject(seu.obj, split.by = "orig.ident"))
seu.list <- lapply(seu.list, DietSeurat, layers = "counts", assays = "RNA")

message(paste0(Sys.time(), " INFO: merging data from k9 and human."))
samNames <- names(seu.list)
seu.merge <- merge(seu.list[1][[1]], y = seu.list[2:length(seu.list)],
                  add.cell.ids = samNames, 
                  project = "hu_k9_comp"
                 )
rm(seu.list)
gc()

#integrate the data
message(paste0(Sys.time(), " INFO: integrating data from k9 and human."))
seu.obj <- integrateData(din = NULL, pattern = NULL,
                          saveRDS = F, 
                          outName = outName,  dout = "../output/s2/",
                          orig.reduction = "pca",
                          method = "HarmonyIntegration", 
                          normalization.method = "LogNormalize", 
                          indReClus = TRUE, seu.obj = seu.merge,
                          runAllMethods = TRUE
                        )
gc()

#complete data visualization & save the RDS file
message(paste0(Sys.time(), " INFO: data integration complete. compeleting dimension reduction and saving integrated object as a .rds file in ../s3/."))
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = outName, 
                        final.dims = 45, final.res = 0.8, stashID = "clusterID", algorithm = 3, min.dist = 0.2, n.neighbors = 20,
                        prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated.harmony",
                        saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "HLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
)
gc()

#update user
message(paste0(Sys.time(), " INFO: file saved."))

