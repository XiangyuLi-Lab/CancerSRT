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

# read data, perform dimension reduction and annotation, just as usual analysis previously
CRC_meta1 <- Load10X_Spatial(data.dir = 'D:\\CRC_GSE206552\\meta1',
                             filename = '/filtered_feature_bc_matrix.h5',
                             assay = "Spatial",
                             slice = 'spatial',
                             filter.matrix = T)
CRC_meta1[["percent.mt"]] <- PercentageFeatureSet(CRC_meta1, pattern = "^MT[-]")
VlnPlot(CRC_meta1,features = c("nCount_Spatial","nFeature_Spatial","percent.mt"))
dim(CRC_meta1@assays$Spatial@counts)
CRC_meta1 <- SCTransform(CRC_meta1, assay="Spatial",verbose = FALSE)
# PCA, tSNE, umap
CRC_meta1 <- RunPCA(CRC_meta1, assay = "SCT", verbose = FALSE)
ElbowPlot(CRC_meta1)
CRC_meta1 <- FindNeighbors(CRC_meta1, reduction = "pca", dims = 1:15)
CRC_meta1 <- FindClusters(CRC_meta1, resolution=0.6,verbose = FALSE)
CRC_meta1 <- RunTSNE(CRC_meta1,reduction = "pca", dims = 1:15)
plot1 <- DimPlot(CRC_meta1, reduction = "tsne", label = T)
CRC_meta1 <- RunUMAP(CRC_meta1, reduction = "pca", dims = 1:15)
plot2 <- DimPlot(CRC_meta1, reduction = "umap", label = T)
# reference
load("D:\\CRC_sc2.RData")
counts <- CRC_sc2@assays$RNA@counts
counts_ref <- CRC_sc2@assays$RNA@counts
colnames(CRC_sc2@meta.data)
pdata <- CRC_sc2@meta.data[,c("celltype","malignancy")]
rownames(pdata) = rownames(CRC_sc2@meta.data)
pdata$Index = NULL
colnames(pdata) <- c("celltype","malignancy")
CRC.ref <- SummarizedExperiment(assays=list(counts=counts),colData = pdata)
CRC.ref <- logNormCounts(CRC.ref)
counts <- CRC_meta1@assays$Spatial@counts
common_gene <- intersect(rownames(counts),rownames(counts_ref))
counts <- counts[common_gene,]
counts_ref <- counts_ref[common_gene,]
CRC_meta1.find <- SummarizedExperiment(assays=list(counts=counts))
CRC_meta1.find <- logNormCounts(CRC_meta1.find)
CRC_meta1.find
CRC_meta1.find1 <- SingleR(test = CRC_meta1.find,ref = CRC.ref,
                           labels = CRC.ref$celltype)
CRC_meta1$celltype <- CRC_meta1.find1$labels
CRC_meta1.find2 <- SingleR(test = CRC_meta1.find,ref = CRC.ref,
                           labels = CRC.ref$malignancy)
CRC_meta1$malignancy <- CRC_meta1.find2$labels

# cellchatv2
color.use <- scPalette(nlevels(CRC_meta1))
names(color.use) <- levels(CRC_meta1)
SpatialDimPlot(CRC_meta1, label = T, label.size = 3, cols = color.use)
data.input = GetAssayData(CRC_meta1, slot = "data", assay = "SCT") # normalized data matrix
meta = data.frame(labels = CRC_meta1@meta.data$celltype, samples = "sample1", row.names = names(Idents(CRC_meta1)))
meta$samples <- factor(meta$samples)
unique(meta$labels)
unique(meta$samples)
spatial.locs = GetTissueCoordinates(CRC_meta1, scale = NULL, cols = c("imagerow", "imagecol"))
scalefactors = jsonlite::fromJSON(txt = file.path("D:\\meta1\\spatial", 'scalefactors_json.json'))
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
min(d.spatial[d.spatial!=0]) 
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
cellchat
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
# if get wrong, check your data if there are 0 or position is not match 
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.range = 250, scale.distance = 0.01,
                              contact.dependent = TRUE, contact.range = 100)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")

# save the information about ligand-receptor for visualization, can be omitted
dim3 <- dim(cellchat@net$prob)[3]
all_zero_planes <- c()
for (i in 1:dim3) {
  if (all(cellchat@net$prob[,,i] == 0)) {
    all_zero_planes <- c(all_zero_planes, i)
  }
}
if (length(all_zero_planes) > 0) {
  cellchat@net$prob <- cellchat@net$prob[,, -all_zero_planes]
}
mdimnames <- dimnames(cellchat@net$prob)
# rowname is ligand
CellL_name <- mdimnames[1][1]
interaction_name <- mdimnames[3][1]
CellR_name <- mdimnames[2][1]
dim3 <- dim(cellchat@net$prob)[3]
dim2 <- dim(cellchat@net$prob)[2]
dim1 <- dim(cellchat@net$prob)[1]
ligand_name <- list()
receptor_name <- list()
pathway_name <- list()
for(i in 1:dim3){
  indice <- grep(interaction_name[[1]][i],cellchat@LR$LRsig$interaction_name)
  ligand_name <- append(ligand_name, list(cellchat@LR$LRsig$ligand[indice]))
  receptor_name <- append(receptor_name, list(cellchat@LR$LRsig$receptor[indice]))
  pathway_name <- append(pathway_name, list(cellchat@LR$LRsig$pathway_name[indice]))
}
ligand_name[[3]]
receptor_name[[5]]
pathway_name[[4]]
mprob <- cellchat@net$prob
res <- data.frame(cellF = character(0), interaction_name = character(0), ligand = character(0), pathway = character(0), receptor = character(0),cellT = character(0))
for(i in 1:dim3){
  for(j in 1:dim2){
    for(k in 1:dim1){
      if(mprob[k,j,i]>0){
        new_row <- data.frame(cellF= CellL_name[[1]][k], interaction_name = interaction_name[[1]][i], ligand = ligand_name[[i]], pathway = pathway_name[[i]], receptor = receptor_name[[i]], cellT = CellR_name[[1]][j])
        res_temp <- rbind(res, new_row)
        res <- res_temp
      }
    }
  }
}
dim(res)
write.csv(res, "D:\\CRC_meta1_cellchat.csv")