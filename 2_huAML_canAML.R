#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/CSUClinHeme/users/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")

#load in processed data from https://doi.org/10.1158/2159-8290.CD-22-1297
seu.obj <- readRDS("../external_data/ssp_mlsc_final.Rds")
seu.obj[["RNA"]] <- as(object = seu.obj[["RNA"]], Class = "Assay5")
seu.obj$orig.ident <- seu.obj$sample
outName <- "human_k9_comp"

#inspect data for proper import -- looks good
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "scArches_Cluster",
              pt.size = 0.1,
              label = T,
              label.box = T,
              repel = F,
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_scArches_UMAP.png"), width = 7, height = 7)


#plot m-LSC score
p <- FeaturePlot(seu.obj, features = "m-LSC score", reduction = "umap") 
p <- formatUMAP(p) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_mLSC_score_UMAP.png"), width = 7, height = 7)


#plot m-LSC score by cluster
VlnPlot(seu.obj, features = c("m-LSC score","p-LSC score"), group.by = "scArches_Cluster", stack = T, flip = T, fill.by = "ident")
ggsave(paste0("../output/", outName, "/", outName, "_mLSC_score_viln.png"), width = 7, height = 5)


### Integrate cross-species

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