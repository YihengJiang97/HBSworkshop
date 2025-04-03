### Before everything 
# https://drive.google.com/file/d/1dzdqGE9Z_K0rcD2IlcMDddRb75SA_z-b/view?usp=drive_link
# install.packages(c("qs", "Seurat", "SCpubr"))

library(Seurat)
library(qs)
library(SCpubr)

## set seed for reproducibility
set.seed(250403)

## set working directory
workdir <- file.path(Sys.getenv("USERPROFILE"), "Desktop") # windows
# wordir <- file.path(Sys.getenv("HOME"), "Desktop") # Mac
setwd(workdir)
getwd()

## create a file to store data
dir.create("workshop0403")
setwd("workshop0403")
getwd()

###################################################################################################
### Step0. ReadData 
###################################################################################################
load("seu_obj_workshop.Rdata")
seurat_obj <- seu_raw

###################################################################################################
### Step1. QC
###################################################################################################
# Add percentage of mitochondrial gene expression
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

## Step1.1. Filter cells
minGene = 200
maxGene = 5000
minUMI = 500
maxUMI = 25000
pctMT = 20


seurat_obj_flter <- subset(seurat_obj, 
                           subset = nFeature_RNA > minGene & 
                             nFeature_RNA < maxGene &
                             nCount_RNA > minUMI &
                             nCount_RNA < maxUMI &
                             percent.mt < pctMT 
)
seurat_obj_flter
rm(seurat_obj, seu_raw)

## Step1.2. Filter genes
count_raw = GetAssayData(object = seurat_obj_flter, assay = "RNA", layer = "counts")
meta_data_raw = seurat_obj_flter@meta.data

seu_obj_new <- CreateSeuratObject(counts = count_raw, min.cells = 3, min.features = 300)
seu_obj_new <- AddMetaData(seu_obj_new, metadata = meta_data_raw)
seu_obj_new

rm(count_raw, meta_data_raw)

###################################################################################################
### Step2. Normalize and Batch Correction
###################################################################################################
seu_obj_new[["RNA"]] <- split(seu_obj_new[["RNA"]], f = seu_obj_new$Sample_id)
seu_obj_new <- NormalizeData(seu_obj_new, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj_new <- FindVariableFeatures(seu_obj_new, selection.method = "vst", nfeatures = 2000)
seu_obj_new <- ScaleData(seu_obj_new, vars.to.regress = c("nCount_RNA", "percent.mt"))
seu_obj_new <- RunPCA(seu_obj_new, features = VariableFeatures(object = seu_obj_new))

### Step2.1. Umap Without Batch Correction
seurat_obj_raw <- seu_obj_new
seurat_obj_raw <- FindNeighbors(seurat_obj_raw, reduction = "pca", dims = 1:10)
seurat_obj_raw <- FindClusters(seurat_obj_raw, resolution = 0.7)
seurat_obj_raw <- RunUMAP(seurat_obj_raw, dims = 1:10)

## plot
pdf("Step2.1.umap_before_Batch_correct.pdf", width = 8, height = 6) 
DimPlot(seurat_obj_raw, reduction = "umap", group.by = c("Sample_id"))
dev.off()

### step2.2 Umap with Batch Correction (Harmony)
seurat_harmony <- IntegrateLayers(
  object = seu_obj_new, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

seurat_harmony <- FindNeighbors(seurat_harmony, reduction = "harmony", dims = 1:10)
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.7)
seurat_harmony <- RunUMAP(seurat_harmony, dims = 1:10, reduction = "harmony", reduction.name = "umap.harmony")

## plot
pdf("Step2.2.umap_after_Batch_correct.pdf", width = 8, height = 6) 
DimPlot(seurat_harmony, reduction = "umap.harmony", group.by = c("Sample_id"))
dev.off()

###################################################################################################
### Step3.Cell Type Annotation
###################################################################################################
Idents(seurat_harmony) <- seurat_harmony$seurat_clusters

genes <- list("Urothelial cells" = c("EPCAM", "KRT19", "CDH1", "KRT18"),
              "T cells" = c("CD3D", "CD3E", "CD3G", "TRAC"),
              "B cells" =  c("CD79A", "IGHM", "IGHG3", "IGHA2"),
              "Myeloid cells" = c("CD68", "MARCO", "FCGR3A", "LYZ"),
              "NK cells" = c("NCAM1", "NKG7", "GNLY", "KLRD1"),
              "MAST cells" = c("KIT", "MS4A2","GATA2"),
              "Fibroblasts" = c("DCN", "COL1A1", "COL1A2", "THY1" ),
              "Endothelial cells" = c("PECAM1", "CLDN5", "FLT1", "RAMP2")
              
)  
# plot
pdf("Step3.1.markergenes_expr.pdf", width = 12, height = 10) 
do_DotPlot(sample = seurat_harmony, features = genes, dot.scale = 12, legend.length = 50,
           legend.framewidth = 2, font.size = 12)
dev.off()

## Assign Cell Types
seurat_harmony$celltype <- plyr::mapvalues(seurat_harmony$seurat_clusters,
                                       from = 0:18,
                                       to = c(
                                         "0" = "T_NK cells",
                                         "1" = "T_NK cells",
                                         "2" = "Fibroblasts",
                                         "3" = "Urothelial cells",
                                         "4" = "Endothelial cells",
                                         "5" = "Myeloid cells",
                                         "6" = "Urothelial cells",
                                         "7" = "Urothelial cells",
                                         "8" = "T_NK cells",
                                         "9" = "B cells",
                                         "10" = "Urothelial cells",
                                         "11" = "T_NK cells",
                                         "12" = "Fibroblasts",
                                         "13" = "Fibroblasts",
                                         "14" = "T_NK cells",
                                         "15" = "Urothelial cells",
                                         "16" = "Endothelial cells",
                                         "17" = "B cells",
                                         "18" = "Fibroblasts"
                                       ))

pdf("Step3.2.umap_celltypes.pdf", width = 15, height = 8) 
DimPlot(seurat_harmony, reduction = "umap.harmony", group.by = c("celltype", "seurat_clusters"), label = T)
dev.off()

### plot3d
# library(plotly)
# seurat_harmony <- RunUMAP(seurat_harmony, dims = 1:10, reduction = "harmony", reduction.name = "umap.harmony", n.components = 3)
# DimPlot(seurat_harmony, reduction = "umap.harmony", group.by = c("celltype"))
# 
# meta_data_plot <- seurat_harmony@meta.data
# umap_embeddings <- Embeddings(seurat_harmony, reduction = "umap.harmony")
# 
# plot.data_3d <- data.frame(
#   row.names = rownames(meta_data_plot),
#   celltype = meta_data_plot$celltype,
#   UMAP_1 = umap_embeddings[,1],
#   UMAP_2 = umap_embeddings[,2],
#   UMAP_3 = umap_embeddings[,3]
# )
# plot.data_3d$label <- paste(rownames(plot.data_3d))
# 
# colors_tmp<-c("lightseagreen","gray50","darkgreen","red4","red","turquoise4","black","yellow4","royalblue1","lightcyan3","peachpuff3","khaki3","gray20","orange2","royalblue4","yellow3","gray80","darkorchid1","lawngreen","plum2","darkmagenta")
# 
# fig <- plot_ly(data = plot.data_3d, 
#              x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
#              color = ~celltype, 
#              # colors = colors_tmp,
#              type = "scatter3d", 
#              mode = "markers", 
#              marker = list(size = 1, width=8), # controls size of points
#              # text=~label, 
#              hoverinfo="text") 
# 
# fig <- fig %>%
#   layout(legend = list(
#     font = list(size = 15),  
#     itemsizing = "constant", 
#     itemwidth = 50           
#   ))
# fig
