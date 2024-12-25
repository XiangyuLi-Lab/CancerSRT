  library(Seurat)
  library(SingleR)
  library(curl)
  library(devtools)
  library(usethis)
  library(tidyverse)
  library(SummarizedExperiment)
  library(scuttle)
  library(patchwork)
  library(CellChat)
  library(NMF)
  library(devtools)
  library(spacexr)
  library(Matrix)
  library(ggplot2)
  library(SeuratData)
  library(copykat)
  library(scCancer)
  library(DropletUtils)
  library(AnnotationDbi)
  library(dplyr)
  library(IRanges)
  library(S4Vectors)
  library(stats)
  library(clusterProfiler)
  library(topGO)
  library(Rgraphviz)
  library(pathview)
  library(org.Hs.eg.db)#mouse
  library(stringr)
  library(graph)
  library(ggnewscale)
  library(dplyr)
  library(msigdbr)
  library(GSVA)
  library(pheatmap)
  library(gelnet)
  library(Rmisc)
  library(tidyverse)
  library(magrittr)


# Data loading (the following code takes Visium_S01_MAR21 as an example) Select the corresponding reading method according to different data.
Visium_S01_MAR21 <- Load10X_Spatial(data.dir = 'yourpath\\Visium_S01_MAR21',
                              filename = '\\filtered_feature_bc_matrix.h5',
                              assay = "Spatial",
                              slice = 'spatial',
                              filter.matrix = T) 
# Data preprocessing
Visium_S01_MAR21[["percent.mt"]] <- PercentageFeatureSet(Visium_S01_MAR21, pattern = "^MT[-]") 
VlnPlot(Visium_S01_MAR21,features = c("nCount_Spatial","nFeature_Spatial","percent.mt"))
dim(Visium_S01_MAR21@assays$Spatial@counts) 

# If an error occurs during the SCT processing, it means that the data has null values and needs to be filtered using this code：  
# Visium_S01_MAR21 <- Visium_S01_MAR21[,Visium_S01_MAR21$nCount_RNA > 0]
Visium_S01_MAR21 <- SCTransform(Visium_S01_MAR21, assay="Spatial" ,verbose = FALSE)

Visium_S01_MAR21 <- RunPCA(Visium_S01_MAR21, assay = "SCT", verbose = FALSE) 
ElbowPlot(Visium_S01_MAR21)
Visium_S01_MAR21 <- FindNeighbors(Visium_S01_MAR21, reduction = "pca", dims = 1:15) 
Visium_S01_MAR21 <- FindClusters(Visium_S01_MAR21, resolution=0.6,verbose = FALSE)
Visium_S01_MAR21 <- RunTSNE(Visium_S01_MAR21,reduction = "pca", dims = 1:15) 
plot1 <- DimPlot(Visium_S01_MAR21, reduction = "tsne", label = T)
Visium_S01_MAR21 <- RunUMAP(Visium_S01_MAR21, reduction = "pca", dims = 1:15)
plot2 <- DimPlot(Visium_S01_MAR21, reduction = "umap", label = T)

# single  
# Used to generate a reference set for subsequent use (the code uses the GBM reference set as an example). If you already have a reference set, you do not need to run this step.
# GBM_sc <- Read10X_h5("yourpath\\GBM\\Glioma_GSE131928_10X_expression.h5")
# GBM_sc <- CreateSeuratObject(counts = GBM_sc, project = "GBM",
#                              min.cells = 3,min.features = 200)
# dim(GBM_sc@assays$RNA@counts)
# GBM_sc[["percent.mt"]] <- PercentageFeatureSet(GBM_sc, pattern = "^MT-")
# VlnPlot(GBM_sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# GBM_sc <- SCTransform(GBM_sc)
# GBM_sc <- RunPCA(GBM_sc)
# ElbowPlot(GBM_sc)
# GBM_sc <- FindNeighbors(GBM_sc, dims = 1:30)
# GBM_sc <- FindClusters(GBM_sc, resolution = 1)
# GBM_sc <- RunUMAP(GBM_sc,dims = 1:30)
# DimPlot(GBM_sc, reduction = "umap")
# GBM_sc <- RunTSNE(GBM_sc,dims = 1:30)
# DimPlot(GBM_sc, reduction = "tsne")
# 
# CellMetainfo_table <- read.delim("yourpath/Glioma_GSE131928_10X_CellMetainfo_table.tsv", row.names=1)
# GBM_sc$celltype <- CellMetainfo_table$Celltype..major.lineage.
# GBM_sc$malignancy <- CellMetainfo_table$Celltype..malignancy.
# table(GBM_sc$celltype)
# table(GBM_sc$malignancy)
# 
# save(GBM_sc,file = "yourpath\\GBM_sc.RData")

# Custom reference set OPSCC_sc.RData ，A reference set generated or downloaded by you
load("yourpath\\OPSCC_sc.RData")
counts <- OPSCC_sc@assays$RNA@counts
counts_ref <- OPSCC_sc@assays$RNA@counts
colnames(OPSCC_sc@meta.data) 
pdata <- OPSCC_sc@meta.data[,c("celltype","malignancy")]  
rownames(pdata) = rownames(OPSCC_sc@meta.data) 
pdata$Index = NULL  

colnames(pdata) <- c("celltype","malignancy") 
OPSCC.ref <- SummarizedExperiment(assays=list(counts=counts),colData = pdata)  
OPSCC.ref <- logNormCounts(OPSCC.ref) 
## SingleR
counts <- Visium_S01_MAR21@assays$Spatial@counts
common_gene <- intersect(rownames(counts),rownames(counts_ref))
# There may be an error here because the gene names in the reference set and the gene names in the expression matrix are inconsistent. 
# The error is reflected in the fact that common_gene in the subsequent code is empty.
# solution：
# This data is provided in the repository for unifying gene names
#data <- read.csv("yourpath\\ensemble.csv")
# Modify the id of all row names
#count = 1  Record execution times
# Record the number of data for which gene names were not found
#no = 1  
#for (i in rownames(counts)){
#   rownames(pdata)[count]<-c(substring(i, 5,))
#   Find other data in the same row based on a specific value
#   row_data <- data[data$ENSEMBL == i, ]
#   other_columns <- row_data[, c("SYMBOL")]
#   if (length(other_columns) == 0) {
#     count = count + 1
#     no = no + 1
#     next
#   }
#   rownames(counts)[count] <- other_columns
#   print(rownames(counts)[count])
#   count = count + 1
#   print(count)
#}
counts <- counts[common_gene,]
counts_ref <- counts_ref[common_gene,]
Visium_S01_MAR21.find <- SummarizedExperiment(assays=list(counts=counts))  
Visium_S01_MAR21.find <- logNormCounts(Visium_S01_MAR21.find)
Visium_S01_MAR21.find  
# Annotate cells
Visium_S01_MAR21.find1 <- SingleR(test = Visium_S01_MAR21.find,ref = OPSCC.ref,   
                            labels = OPSCC.ref$celltype)
Visium_S01_MAR21$celltype <- Visium_S01_MAR21.find1$labels 

Visium_S01_MAR21.find2 <- SingleR(test = Visium_S01_MAR21.find,ref = OPSCC.ref, 
                            labels = OPSCC.ref$malignancy)
Visium_S01_MAR21$malignancy <- Visium_S01_MAR21.find2$labels

## cellchat 
# Conduct cell-cell communication analysis
table(Visium_S01_MAR21$celltype)
table(Visium_S01_MAR21$malignancy)
cellchat <- createCellChat(Visium_S01_MAR21,group.by = "celltype") 
groupsize <- as.numeric(table(Visium_S01_MAR21$celltype)) 
# Using the Human Basic Database
cellchatDB <- CellChatDB.human
use <- subsetDB(cellchatDB,search = "Secreted Signaling") 
cellchat@DB <- use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat) 
cellchat <- computeCommunProb(cellchat,raw.use = TRUE,population.size = TRUE) 
cellchat <- filterCommunication(cellchat,min.cells = 0) 
df.net <-  subsetCommunication(cellchat) 
cellchat <- computeCommunProbPathway(cellchat) 
df.netp <- subsetCommunication(cellchat,slot.name = "netP")  
cellchat <- aggregateNet(cellchat) 
groupSize <- as.numeric(table(cellchat@idents)) 
# Ensure that the output results are always below this result set
dir.create("yourpath/Visium_S01_MAR21")
setwd("yourpath/Visium_S01_MAR21")  
netVisual_circle(cellchat@net$count,vertex.weight = groupSize,weight.scale = T,
                 arrow.width = 0.5,arrow.size = 0.5,
                 label.edge = F,title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight,vertex.weight = groupSize,weight.scale = T,
                 arrow.width = 0.5,arrow.size = 0.5,
                 label.edge = F,title.name = "Interaction weights/strength")
## RCTD 
# Deconvolution analysis was used to identify cell types in the input matrix.
sc <- OPSCC_sc
counts <- as.matrix(sc@assays$RNA@counts)
counts <- apply(counts, 2, as.integer)
rownames(counts) <- rownames(sc@assays$RNA@counts)
meta_data <- sc@meta.data
cell_types <- meta_data$celltype; names(cell_types) <- rownames(meta_data)
cell_types <- as.factor(cell_types)
nUMI <- meta_data$nCount_RNA
nUMI <- as.integer(nUMI)
names(nUMI) <- rownames(meta_data)
cell_types <- gsub("/",".",cell_types)
cell_types <- as.factor(cell_types)
names(cell_types) <- rownames(meta_data)

### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)
table(reference@cell_types)
## 空间转录组
spatial <- Visium_S01_MAR21
spatialcounts <- as.matrix(spatial@assays$Spatial@counts)
spatialcounts <- apply(spatialcounts, 2, as.integer)
rownames(spatialcounts) <- rownames(spatial@assays$Spatial@counts)
coords <- spatial@images$spatial@coordinates[,2:3]
nUMI <- spatial@meta.data$nCount_SCT
names(nUMI)<- rownames(spatial@meta.data)
rownames(coords) <- colnames(spatialcounts) <- names(nUMI)
### Create SpatialRNA object
puck <- SpatialRNA(coords, spatialcounts, nUMI)

# Result prediction
myRCTD <- create.RCTD(puck, reference,max_cores = 8,CELL_MIN_INSTANCE = 0)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet') 


# Results Visualization
results <- myRCTD@results
results
norm_weights <- sweep(results$weights,1,rowSums(results$weights),"/")
RCTD_mat <- data.frame(norm_weights)
write.csv(RCTD_mat,"yourpatht\\Visium_S01_MAR21\\Visium_S01_MAR21_RCTD.csv")

##infercnv
dir.create("yourpath\\Visium_S01_MAR21\\scCancer")
exp <- Visium_S01_MAR21@assays$Spatial@counts
gene <- rownames(Visium_S01_MAR21@assays$Spatial@counts)
EntrezID <- bitr(gene, fromType="SYMBOL", toType="ENSEMBL", OrgDb ="org.Hs.eg.db",drop = TRUE)
index <- EntrezID[!duplicated(EntrezID$SYMBOL),]
rownames(index) <- index$SYMBOL
index <- index[order(index$SYMBOL),]
index <- na.omit(index)
exp <- exp[rownames(index),]
index$gene <- index$SYMBOL
dim(index)
# A path to save the results
generate10Xdata(exp,index[,2:3],outPath = "yourpath\\Visium_S01_MAR21\\scCancer")
dataPath <- "yourpath\\Visium_S01_MAR21\\scCancer"     
statPath <- "yourpath\\Visium_S01_MAR21\\scCancer" 
savePath <- "yourpath\\Visium_S01_MAR21\\scCancer"  
# The sample name
sampleName <- "Visium_S01_MAR21"          
# The author name to mark the report
authorName <- "huoyuying"           
# Generate report-figures cellQC cekkMan geneMan report-sc
stat.results <- runScStatistics(
  dataPath = dataPath,
  savePath = savePath,
  authorName = authorName,
  sampleName = sampleName,
) 
# Generate annotated cluster expr rds file
anno.results <- runScAnnotation(
  dataPath = dataPath,
  statPath = statPath,
  savePath = savePath,
  authorName = authorName,
  sampleName = sampleName
)

expr <- readRDS("yourpath\\Visium_S01_MAR21\\scCancer\\expr.RDS")
pred <- data.frame(expr@meta.data$Malign.type)
table(pred$expr.meta.data.Malign.type)
rownames(pred) <- rownames(expr@meta.data)
colnames(pred) <- "label"
name <- gsub("-1","",rownames(Visium_S01_MAR21@meta.data))
pred <- pred[name,]
Visium_S01_MAR21$inferCNV.pred <- pred
table(Visium_S01_MAR21$inferCNV.pred)

#copykat
dir.create("yourpath\\Visium_S01_MAR21\\copykat")
setwd("yourpath\\Visium_S01_MAR21\\copykat")
exp <- as.matrix(Visium_S01_MAR21@assays$Spatial@counts)
copykat <- copykat(exp,id.type = "S", 
                   ngene.chr = 5,  
                   win.size = 25, 
                   # Increasing KS.cut will reduce sensitivity, usually in the range of 0.05-0.15
                   KS.cut = 0.1,
                   sam.name = "test",distance = "euclidean", norm.cell.names = "",n.cores = 1)

pred <- data.frame(copykat$prediction)
CNA <- data.frame(copykat$CNAmat)

table(pred$copykat.pred)
pred$copykat.pred <- gsub("c1:diploid:low.conf","diploid",pred$copykat.pred)
pred$copykat.pred <- gsub("c2:aneuploid:low.conf","aneuploid",pred$copykat.pred)
pred$copykat.pred <- gsub("not.defined","unknown",pred$copykat.pred)

pred <- pred[rownames(Visium_S01_MAR21@meta.data),]
Visium_S01_MAR21$copykat.pred <- pred$copykat.pred
table(Visium_S01_MAR21$copykat.pred)
# 14 functional status analysis 
# The signatures.csv file is provided in the repository
signatures <- read.csv("yourpath\\signatures.csv")
expr <- as.data.frame(Visium_S01_MAR21@assays$Spatial@data) 
expr = as.matrix(expr)
msigdbr_list = split(x = signatures$GeneName, f = signatures$Signature)
res <- gsva(expr, msigdbr_list, kcdf="Gaussian",method = "gsva",parallel.sz=10) 
# Angiogenesis
Visium_S01_MAR21$Angiogenesis <- res[1,]  
temp <- (Visium_S01_MAR21$Angiogenesis-min(Visium_S01_MAR21$Angiogenesis))/(max(Visium_S01_MAR21$Angiogenesis)-min(Visium_S01_MAR21$Angiogenesis))
Visium_S01_MAR21$Angiogenesis <- temp
# Apoptosis
Visium_S01_MAR21$Apoptosis <- res[2,]
temp <- (Visium_S01_MAR21$Apoptosis-min(Visium_S01_MAR21$Apoptosis))/(max(Visium_S01_MAR21$Apoptosis)-min(Visium_S01_MAR21$Apoptosis))
Visium_S01_MAR21$Apoptosis <- temp
# CellCycle
Visium_S01_MAR21$CellCycle <- res[3,]
temp <- (Visium_S01_MAR21$CellCycle-min(Visium_S01_MAR21$CellCycle))/(max(Visium_S01_MAR21$CellCycle)-min(Visium_S01_MAR21$CellCycle))
Visium_S01_MAR21$CellCycle <- temp
# Differentiation
Visium_S01_MAR21$Differentiation <- res[4,]
temp <- (Visium_S01_MAR21$Differentiation-min(Visium_S01_MAR21$Differentiation))/(max(Visium_S01_MAR21$Differentiation)-min(Visium_S01_MAR21$Differentiation))
Visium_S01_MAR21$Differentiation <- temp
# DNADamage
Visium_S01_MAR21$DNADamage <- res[5,]
temp <- (Visium_S01_MAR21$DNADamage-min(Visium_S01_MAR21$DNADamage))/(max(Visium_S01_MAR21$DNADamage)-min(Visium_S01_MAR21$DNADamage))
Visium_S01_MAR21$DNADamage <- temp
# DNARepair
Visium_S01_MAR21$DNARepair <- res[6,]
temp <- (Visium_S01_MAR21$DNARepair-min(Visium_S01_MAR21$DNARepair))/(max(Visium_S01_MAR21$DNARepair)-min(Visium_S01_MAR21$DNARepair))
Visium_S01_MAR21$DNARepair <- temp
# EMT
Visium_S01_MAR21$EMT <- res[7,]
temp <- (Visium_S01_MAR21$EMT-min(Visium_S01_MAR21$EMT))/(max(Visium_S01_MAR21$EMT)-min(Visium_S01_MAR21$EMT))
Visium_S01_MAR21$EMT <- temp
# Hypoxia
Visium_S01_MAR21$Hypoxia <- res[8,]
temp <- (Visium_S01_MAR21$Hypoxia-min(Visium_S01_MAR21$Hypoxia))/(max(Visium_S01_MAR21$Hypoxia)-min(Visium_S01_MAR21$Hypoxia))
Visium_S01_MAR21$Hypoxia <- temp
# Inflammation
Visium_S01_MAR21$Inflammation <- res[9,]
temp <- (Visium_S01_MAR21$Inflammation-min(Visium_S01_MAR21$Inflammation))/(max(Visium_S01_MAR21$Inflammation)-min(Visium_S01_MAR21$Inflammation))
Visium_S01_MAR21$Inflammation <- temp
# Invasion
Visium_S01_MAR21$Invasion <- res[10,]
temp <- (Visium_S01_MAR21$Invasion-min(Visium_S01_MAR21$Invasion))/(max(Visium_S01_MAR21$Invasion)-min(Visium_S01_MAR21$Invasion))
Visium_S01_MAR21$Invasion <- temp
# Metastasis
Visium_S01_MAR21$Metastasis <- res[11,]
temp <- (Visium_S01_MAR21$Metastasis-min(Visium_S01_MAR21$Metastasis))/(max(Visium_S01_MAR21$Metastasis)-min(Visium_S01_MAR21$Metastasis))
Visium_S01_MAR21$Metastasis <- temp
# Proliferation
Visium_S01_MAR21$Proliferation <- res[12,]
temp <- (Visium_S01_MAR21$Proliferation-min(Visium_S01_MAR21$Proliferation))/(max(Visium_S01_MAR21$Proliferation)-min(Visium_S01_MAR21$Proliferation))
Visium_S01_MAR21$Proliferation <- temp
# Quiescence
Visium_S01_MAR21$Quiescence <- res[13,]
temp <- (Visium_S01_MAR21$Quiescence-min(Visium_S01_MAR21$Quiescence))/(max(Visium_S01_MAR21$Quiescence)-min(Visium_S01_MAR21$Quiescence))
Visium_S01_MAR21$Quiescence <- temp
# Stemness
Visium_S01_MAR21$Stemness <- res[14,]
temp <- (Visium_S01_MAR21$Stemness-min(Visium_S01_MAR21$Stemness))/(max(Visium_S01_MAR21$Stemness)-min(Visium_S01_MAR21$Stemness))
Visium_S01_MAR21$Stemness <- temp

# Output File
Visium_S01_MAR21_coord <- Visium_S01_MAR21@images$spatial@coordinates
Visium_S01_MAR21_umap <- data.frame(FetchData(Visium_S01_MAR21, vars=c("UMAP_1", "UMAP_2")))
Visium_S01_MAR21_tsne <- data.frame(FetchData(Visium_S01_MAR21, vars=c("tSNE_1", "tSNE_2")))
Visium_S01_MAR21_anno <- Visium_S01_MAR21@meta.data
# The above analysis may not be performed. Infercnv and copycat may not be performed during the analysis. 
# Therefore, you need to delete some content in the output.
Visium_S01_MAR21_anno <- Visium_S01_MAR21@meta.data[,c("seurat_clusters","celltype","malignancy",
                                           "inferCNV.pred","copykat.pred",
                                           "Angiogenesis","Apoptosis","CellCycle","Differentiation",
                                           "DNADamage","DNARepair","EMT","Hypoxia","Inflammation",
                                           "Invasion","Metastasis","Proliferation",
                                           "Quiescence","Stemness")]
Visium_S01_MAR21_anno <- cbind(Visium_S01_MAR21_coord,Visium_S01_MAR21_umap,Visium_S01_MAR21_tsne,Visium_S01_MAR21_anno)
write.csv(Visium_S01_MAR21_anno,"yourpath\\Visium_S01_MAR21\\Visium_S01_MAR21_anno.csv")

## markers
Idents(Visium_S01_MAR21) <- Visium_S01_MAR21@meta.data$seurat_clusters
Visium_S01_MAR21$orig.ident <- Visium_S01_MAR21@meta.data$seurat_clusters
data.markers <- FindAllMarkers(Visium_S01_MAR21,only.pos = FALSE,min.pct = 0.25,logfc.threshold = 0.25,test.use = "wilcox")
write.csv(data.markers,"yourpath\\Visium_S01_MAR21\\Visium_S01_MAR21_markers_domain.csv",row.names =F)

Idents(Visium_S01_MAR21) <- Visium_S01_MAR21@meta.data$celltype
Visium_S01_MAR21$orig.ident <- Visium_S01_MAR21@meta.data$celltype
data.markers <- FindAllMarkers(Visium_S01_MAR21,only.pos = FALSE,min.pct = 0.25,logfc.threshold = 0.25,test.use = "wilcox")
write.csv(data.markers,"yourpath\\Visium_S01_MAR21\\Visium_S01_MAR21_markers_celltype.csv",row.names =F)

save(Visium_S01_MAR21,Visium_S01_MAR21_anno,RCTD_mat,
     file = "yourpath\\Visium_S01_MAR21\\Visium_S01_MAR21.RData")

## spatialDE
sp_coord <- Visium_S01_MAR21@images$spatial@coordinates[,2:3]
dim(sp_coord)
dim(Visium_S01_MAR21@assays[["Spatial"]]@counts)
## coord
cell <- rownames(sp_coord)
x <- sp_coord$row
y <- sp_coord$col
coord <- data.frame(x,y)
rownames(coord)  <- cell
##counts
counts <- data.frame(Visium_S01_MAR21@assays[["Spatial"]]@counts)
total_counts <- data.frame(colSums(counts))

coord <- data.frame(coord,total_counts)
colnames(coord) <- c("x","y","total_counts")

colnames(counts) <- rownames(coord)

dir.create("yourpath\\Visium_S01_MAR21\\SV")
setwd("yourpath\\Visium_S01_MAR21\\SV")
write.csv(coord,"Visium_S01_MAR21_coord_sp.csv")
write.csv(counts,"Visium_S01_MAR21_counts_sp.csv")

# Immunoassay
immune <- read.csv("yourpath\\immune.csv")
table(immune)
celltype = c("B","CD4Tconv","CD8T","CD8Tex","DC",
             "Endothelial","Fibroblasts","Mast","Mono/Macro","Myofibroblasts",
             "Neutrophils","NK","pDC","Plasma","TMKI67","Treg")
for (i in 1:length(celltype)) {
  gene_set <- immune[immune$celltype == celltype[i],]
  Visium_S01_MAR21 <- AddModuleScore(Visium_S01_MAR21, 
                               features = list(gene_set$markers),
                               name = celltype[i])
  print(celltype[i])
}
colnames(Visium_S01_MAR21@meta.data)
# The cell type may change during the analysis. If an error message appears, the parameters need to be modified.
Visium_S01_MAR21_immune <- Visium_S01_MAR21@meta.data[,c("B1","CD4Tconv1","CD8T1","CD8Tex1","DC1",
                                             "Endothelial1","Fibroblasts1","Mast1","Mono.Macro1","Myofibroblasts1",
                                             "Neutrophils1","NK1","pDC1","Plasma1","TMKI671","Treg1")]
colnames(Visium_S01_MAR21_immune) <- celltype
Visium_S01_MAR21_immune <- cbind (Visium_S01_MAR21_coord,Visium_S01_MAR21_umap,Visium_S01_MAR21_tsne,Visium_S01_MAR21_immune)
write.csv(Visium_S01_MAR21_immune,"yourpath\\Visium_S01_MAR21\\Visium_S01_MAR21_immune.csv")

# TLS Analysis
obj = Visium_S01_MAR21
markers = c("FDCSP","CR2","CXCL13","LTF","CD52","MS4A1","CCL19","LINC00926","LTB","CORO1A",
            "CD79B","TXNIP","CD19","LIMD2","CD37","ARHGAP45","BLK","TMC8","CCL21","PTPN6",
            "ATP2A3","IGHM","SPIB","TMSB4X","CXCR4","NCF1","CD79A","ARHGAP9","DEF6","EVL",
            "TBC1D10C","RASAL3","INPP5D","RNASET2","RASGRP2","TNFRSF13C","RAC2","CD22","ARHGEF1","AC103591.3",
            "TRAF3IP3","HLA-DQB1","CD53","ARHGAP4","TRBC2","POU2AF1","TRAF5","OGA","FCRL3","HLA-DQA1")
temp <- GetAssayData(obj, slot = "data")
temp <- as.data.frame(temp)
temp <- temp[markers,]
temp[is.na(temp)] <- 0
Visium_S01_MAR21$TLS_score <- Matrix::colMeans(temp)
Visium_S01_MAR21$TSL_tag <- NA
max_value <- max(Visium_S01_MAR21$TLS_score)
Visium_S01_MAR21$TSL_tag <- ifelse(Visium_S01_MAR21$TLS_score >= 0.6 * max_value, 'Yes', 'No')
table(Visium_S01_MAR21$TSL_tag)
SpatialFeaturePlot(Visium_S01_MAR21, features = "TLS_score",pt.size.factor = 2)
exp_marker <- Visium_S01_MAR21@meta.data[,c("TLS_score","TSL_tag")]
write.csv(exp_marker,"yourpath\\Visium_S01_MAR21\\Visium_S01_MAR21_TLS.csv")
