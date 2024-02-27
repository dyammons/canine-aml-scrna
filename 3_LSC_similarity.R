#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/CSUClinHeme/users/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")

#load in the canine data
seu.obj.k9 <- readRDS("../output/s3/can_aml_res0.6_dims45_dist0.1_neigh10_S3.rds")

#load in the human data - reprocessed
seu.obj.hu <- readRDS("../output/s3/ssp_mlsc_res0.6_dims45_dist0.1_neigh10_S3.rds")

#set output paramaters
outName <- human_k9_comp

### Example code to get you started (probs dont't wanna change above this line) ###

#plot b7-h3
features <- c("CD276")
p1 <- prettyFeats(seu.obj = seu.obj.hu, nrow = 1, ncol = 1, title.size = 14, features = features, 
                  order = T, legJust = "top", reduction = "umap.integrated", noLegend = T) 
ggsave(paste0("../output/", outName, "/", outName, "_hu_b7h3_featPlots.png"), width = 4, height = 4)

p2 <- prettyFeats(seu.obj = seu.obj.k9, nrow = 1, ncol = 1, title.size = 14, features = features, 
                  order = T, legJust = "top", reduction = "umap.integrated", noLegend = T)
ggsave(paste0("../output/", outName, "/", outName, "_can_b7h3_featPlots.png"), width = 4, height = 4)


#calculate defining features of each canine cluster
vilnPlots(seu.obj = seu.obj.k9, inFile = NULL, groupBy = "celltype.l2", numOfFeats = 24, outName = "canine",
                      outDir = "../output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T, resume = F, resumeFile = NULL
                     )

#calculate defining features of each human cluster
seu.obj.hu$celltype.l2_merged <- droplevels(seu.obj.hu$celltype.l2_merged)
vilnPlots(seu.obj = seu.obj.hu, inFile = NULL, groupBy = "celltype.l2_merged", numOfFeats = 24, outName = "human",
                      outDir = "../output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T, resume = F, resumeFile = NULL, returnViln = F
                     )

#load in cell type gene signatures
dog.df <- read.csv("../output/viln/canine_gene_list.csv")
human.df <- read.csv("../output/viln/human_gene_list.csv")

#check the data quality
dog.df %>% group_by(cluster) %>% summarize(nn = n()) %>% summarize(min = min(nn))
human.df %>% group_by(cluster) %>% summarize(nn = n()) %>% summarize(min = min(nn))

#calc the jaccard similarity index
res <- lapply(unique(dog.df$cluster), function(x){
    dog.list <- dog.df[dog.df$cluster == x, ] %>% .$gene
    
    res_pre <- lapply(unique(human.df$cluster), function(y){
        human.list <- human.df[human.df$cluster == y, ] %>% .$gene
        
        interSect <- length(intersect(human.list, dog.list)) 
        uni <- length(human.list) + length(dog.list) - interSect 
        JaccardIndex <- interSect/uni
        names(JaccardIndex) <- x
        
        return(JaccardIndex)
    })
    
    res_pre <- do.call(rbind, res_pre)
    rownames(res_pre) <- unique(human.df$cluster)
    return(res_pre)
    
})
    
res1 <- do.call(cbind, res)

# #optionally order the celltypes
# rowTarg <- c("Endothelial","Fibroblast" ,"Osteoblast","Cycling osteoblast",
#        "CD14_monocyte","NR4A3_Macrophage","TXNIP_Macrophage","FABP5_Macrophage",
#        "IFN-TAM","Pre-OC","Mature-OC","pDC","cDC2","Mast",
#        "CD8 T cell","CD4 T cell","IFN-sig T cell","NK cell",
#        "Plasma cell","B cell")
# colTarg <- c("Endothelial cell","Fibroblast","Osteoblast_1","Osteoblast_2","Osteoblast_3","Hypoxic_osteoblast",
#          "IFN-osteoblast","Osteoblast_cycling","CD4+_TIM","CD4-_TIM","ANGIO_TAM","TAM_INT","TAM_ACT","LA-TAM_SPP2_hi",
#          "LA-TAM_C1QC_hi","IFN-TAM","Cycling_OC","CD320_OC","Mature_OC","pDC","cDC1","cDC2","mregDC",
#          "Mast cell","Neutrophil","CD8 T cell","CD4 T cell","T_IFN","T_cycling","NK","Plasma cell","B cell")
# res1 <- res1[match(rowTarg, rownames(res1)),]        
# res1 <- res1[ ,match(colTarg, colnames(res1))]   
       
#plot the data
png(file = paste0("../output/", outName, "/", outName, "_jaccard.png"), width=4000, height=4000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(t(res1), #name = "mat", #col = col_fun,
              name = "Jaccard similarity index",
              cluster_rows = F,
              row_title = "Canine cell types",
              row_title_gp = gpar(fontsize = 24),
              col=viridis(option = "magma",100),
              cluster_columns = F,
              column_title = "Human cell types",
              column_title_gp = gpar(fontsize = 24),
              column_title_side = "bottom",
              column_names_rot = 45,
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
             )
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),heatmap_legend_side = "top")
dev.off()