#!/usr/bin/Rscript

### Analysis note:
# This scripts preprocess canine data (with conversion of gene symbols to human)
# and re-preprocesses the human AML dataset to enable label transfer.

#load custom functions & packages
source("/pl/active/CSUClinHeme/users/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")
outName <- "can_aml"

#read in processed k9 data
seu.obj.k9 <- readRDS("../output/s3/outputcombined_BM_CD34_AML_res0.6_dims50_S3.rds")
cnts <- seu.obj.k9@assays$RNA$counts
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "dog",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))
seu.obj.k9 <- CreateSeuratObject(cnts, project = "humanConvert", assay = "RNA",
                                  min.cells = 0, min.features = 0, names.field = 1,
                                  names.delim = "_", meta.data = seu.obj.k9@meta.data)

#split then merge objects
message(paste0(Sys.time(), " INFO: splitting data from k9 object containing human gene symbols."))
seu.list <- SplitObject(seu.obj.k9, split.by = "orig.ident")
seu.list <- lapply(seu.list, DietSeurat, layers = "counts", assays = "RNA")

message(paste0(Sys.time(), " INFO: merging data."))
samNames <- names(seu.list)
seu.merge <- merge(seu.list[1][[1]], y = seu.list[2:length(seu.list)],
                  add.cell.ids = samNames, 
                  project = outName
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
                          runAllMethods = FALSE
                        )
gc()

#complete data visualization & save the RDS file
message(paste0(Sys.time(), " INFO: data integration complete. compeleting dimension reduction and saving integrated object as a .rds file in ../s3/."))
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = outName, 
                       final.dims = 45, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
                       prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated",
                       saveRDS = T, return_obj = T, returnFeats = T,
                       features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                    "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                    "CD4", "MS4A1", "PPBP","HBM")
                      )

gc()

#prepare the human dataset
outName <- "ssp_mlsc"
hu.reference <- readRDS("../external_data/ssp_mlsc_final.Rds")
hu.reference[["RNA"]] <- as(object = hu.reference[["RNA"]], Class = "Assay5")
hu.reference$orig.ident <- hu.reference$sample

#split then merge objects
message(paste0(Sys.time(), " INFO: splitting data from k9 object containing human gene symbols."))
seu.list <- SplitObject(hu.reference, split.by = "orig.ident")
seu.list <- lapply(seu.list, DietSeurat, layers = "counts", assays = "RNA")

message(paste0(Sys.time(), " INFO: merging data."))
samNames <- names(seu.list)
seu.merge <- merge(seu.list[1][[1]], y = seu.list[2:length(seu.list)],
                  add.cell.ids = samNames, 
                  project = outName
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
                          runAllMethods = FALSE
                        )
gc()

#complete data visualization & save the RDS file
message(paste0(Sys.time(), " INFO: data integration complete. compeleting dimension reduction and saving integrated object as a .rds file in ../s3/."))
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = outName, 
                       final.dims = 45, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.1, n.neighbors = 10,
                       prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated",
                       saveRDS = T, return_obj = T, returnFeats = T,
                       features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                    "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                    "CD4", "MS4A1", "PPBP","HBM")
                      )

gc()


### Transfer labels over from https://doi.org/10.1158/2159-8290.CD-22-1297
seu.obj <- readRDS("../output/s3/can_aml_res0.6_dims45_dist0.1_neigh10_S3.rds")
hu.reference <- readRDS("../output/s3/ssp_mlsc_res0.6_dims45_dist0.1_neigh10_S3.rds")
outName <- "human_k9_comp"

#transfer scArches_Cluster annotations
ref.anchors <- FindTransferAnchors(reference = hu.reference, query = seu.obj, dims = 1:30,
    reference.reduction = "pca", features = rownames(seu.obj)[rownames(seu.obj) %in% rownames(hu.reference)])
predictions <- TransferData(anchorset = ref.anchors, refdata = hu.reference$scArches_Cluster, dims = 1:30)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)


#plot the previously annotated labels (only for healthy cells)
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated", 
              group.by = "minorIdent",
              pt.size = 0.1,
              label = T,
              label.box = T,
              repel = F
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_minorIdent_UMAP.png"), width = 7, height = 7)


#plot transfered labels
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated", 
              group.by = "predicted.id",
              pt.size = 0.1,
              label = T,
              label.box = T,
              repel = F
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_scArches_transfer_UMAP.png"), width = 7, height = 7)


#transfer lsc_class annotations
predictions <- TransferData(anchorset = ref.anchors, refdata = hu.reference$lsc_class, dims = 1:30)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)


seu.obj$predicted.id <- factor(seu.obj$predicted.id)
#inspect data for proper import -- looks good
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated", 
              group.by = "predicted.id",
              split.by = "predicted.id",
              pt.size = 0.1,
              ncol = 3,
              label = F,
              label.box = F,
              repel = F
) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_lsc_class_transfer_UMAP.png"), width = 12, height = 4)

