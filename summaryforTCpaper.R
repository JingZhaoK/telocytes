
#### Combine all 10x datasets ####
L120b=Read10X_h5("L120b_output_filtered.h5")
L120a=Read10X_h5('L120a_output_filtered.h5')
SHAM=Read10X_h5("SHAM-CTRL_output_filtered.h5")
colnames(L120b)=paste0("L120b_",colnames(L120b))
colnames(L120a)=paste0("L120a_",colnames(L120a))
colnames(SHAM)=paste0("SHAM-CTRL_",colnames(SHAM))
cnts=cbind(L120b,L120a,SHAM)
dim(cnts)
s=CreateSeuratObject(cnts)

#### QC for all 10x ####
s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^mt-")
VlnPlot(s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
library(patchwork)
plot1 <- FeatureScatter(s, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
s <- subset(s, subset = nFeature_RNA > 1000 & percent.mt < 20)


#### Dataset integration #####

s@meta.data$batch=sapply(strsplit(colnames(s),"_"),"[[",1)

s.list <- SplitObject(s, split.by = "batch")

s.list <- lapply(X = s.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

s.anchors <- FindIntegrationAnchors(object.list = s.list, dims = 1:20)
s.combined <- IntegrateData(anchorset = s.anchors, dims = 1:20)
DefaultAssay(s.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
s.combined <- ScaleData(s.combined, verbose = FALSE)
s.combined <- RunPCA(s.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
s.combined <- RunUMAP(s.combined, reduction = "pca", dims = 1:30)
s.combined <- FindNeighbors(s.combined, reduction = "pca", dims = 1:30)
s.combined <- FindClusters(s.combined, resolution = 1)

#Names of levels (@batch)
Idents(s.combined) <- s.combined$batch
s.combined <-RenameIdents(s.combined,'L120a'='Ligature','L120b'='Ctrl','SHAM-CTRL'='Ctrl')
s.combined$batch <- s.combined@active.ident
p1 <- DimPlot(s.combined, reduction = "umap", group.by = "batch")
Idents(s.combined) <- s.combined$seurat_clusters
p2 <- DimPlot(s.combined, reduction = "umap", label = TRUE)
library(cowplot)
theme_set(theme_cowplot())
plot_grid(p1, p2)
markers<-FindAllMarkers(s.combined)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(s.combined, features = top10$gene) + NoLegend()

#Annotation
markers %>% filter(cluster==0) %>% top_n(n = 10, wt = avg_logFC)
Idents(s.combined) <- s.combined$seurat_clusters
s.combined <- RenameIdents(s.combined, 
                           '0'= "Endothelial cells",
                           '1' = 'Macrophages 1', 
                           '2'= 'Fibroblasts',
                           '3' ="Macrophages 2",
                           '4' = 'Neutrophils 2',
                           '5' = "Intermediate cells", 
                           '6' = 'Proliferative cells',
                           '7' = 'Cemento/osteoblasts', 
                           '8' = "Macrophages 3",
                           '9' = "Epithelial cells",
                           '10'='B cells',
                           '11'='Dendritic cells',
                           '12'='Macrophages 4',
                           '13'= 'Neutrophils 1',
                           '14'='Telocytes',
                           '15'='Fibroblasts',
                           '16'='T cells',
                           '17'="Neutrophils 3")
s.combined$my_clusters<-Idents(s.combined) 
DimPlot(s.combined)
Idents(s.combined) <- s.combined$my_clusters
DimPlot(s.combined,label = T)+NoLegend()
my_levels <- c("Endothelial cells",
               "Epithelial cells",
               'Proliferative cells',
               'Telocytes',
               'Fibroblasts',
               "Intermediate cells",
               'Cemento/osteoblasts',
               'Neutrophils 1',
               'Neutrophils 2',
               'Neutrophils 3',
               'Macrophages 1',
               "Macrophages 2",
               "Macrophages 3",
               'Macrophages 4',
               'Dendritic cells',
               'T cells',
               'B cells'
)
s.combined@active.ident <- factor (x= s.combined@active.ident,levels=my_levels)
#For heatmap-use seurat cluster
Idents(s.combined) <- s.combined$seurat_clusters
DoHeatmap(s.combined, features = top10$gene) 
p<-DoHeatmap(s.combined, features = top10$gene)
DimPlot(s.combined)

sdata<- subset(s.combined, idents = c('Macrophages 1','Macrophages 2','Macrophages 3','Macrophages 4')) #s is combined dataset

# Run the standard workflow for visualization and clustering
sdata <- ScaleData(sdata, verbose = FALSE)
sdata<- FindVariableFeatures(sdata) 
sdata <- RunPCA(sdata, npcs = 15, verbose = FALSE)
# t-SNE and Clustering
sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:15)
sdata <- FindNeighbors(sdata, reduction = "pca", dims = 1:15)
sdata<- FindClusters(sdata, resolution = 0.5)
library(cowplot)
p1 <- DimPlot(sdata, reduction = "umap", group.by = "batch")
p2 <- DimPlot(sdata, reduction = "umap", label = TRUE)+ggthemes::scale_color_tableau()
plot_grid(p1, p2)


###### Cell Chat #####
install.packages('devtools')
install.packages("BiocManager")
BiocManager::install("Biobase")
devtools::install_github("sqjin/CellChat")
#https://github.com/sqjin/CellChat/blob/master/vignettes/CellChat-vignette.Rmd
#https://www.ershicimi.com/p/5df074c82bf1f00302944b37d649467a
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
options(stringsAsFactors = FALSE)

######CTRL######

data.input.ctrl= Ctrl@assays$RNA@data
identity.ctrl <- data.frame(group=Ctrl$my_clusters,row.names = names(Ctrl$my_clusters))
unique(identity$group)
cellchat.ctrl<- createCellChat(data = data.input.ctrl)
cellchat.ctrl
#Add cell information into meta slot of the object
cellchat.ctrl<-addMeta(cellchat.ctrl,meta = identity.ctrl,meta.name = 'labels')
cellchat.ctrl<-setIdent(cellchat.ctrl,ident.use = 'labels')
levels(cellchat.ctrl@idents)
cellchat.ctrl@idents<-factor(x=cellchat.ctrl@idents,levels=my_levels)
groupsize.ctrl<- as.numeric(table(cellchat.ctrl@idents))


#Set the ligand-receptor interaction database

CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$ißnteraction)
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$complex)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
#Select which database to use:
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # use Secreted Signaling for cell-cell communication analysis
cellchat.ctrl@DB <- CellChatDB.use # set the used database in the object

#Preprocessing the expression data for cell-cell communication analysis
cellchat.ctrl <- subsetData(cellchat.ctrl) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat.ctrl <- identifyOverExpressedGenes(cellchat.ctrl)
cellchat.ctrl <- identifyOverExpressedInteractions(cellchat.ctrl)
cellchat.ctrl <- projectData(cellchat.ctrl, PPI.mouse)

#Inference of cell-cell communication network ##
#Compute the communication probability and infer cellular communication network
cellchat.ctrl <- computeCommunProb(cellchat.ctrl)
#fer the cell-cell communication at a signaling pathway level
cellchat.ctrl <- computeCommunProbPathway(cellchat.ctrl)
#Calculate the aggregated cell-cell communication network
cellchat.ctrl <- aggregateNet(cellchat.ctrl)

######Visualization and systems analysis of cell-cell communication network
#Visualize each signaling pathway using hierarchy plot or circle plot

pathways.show <- c('HGF') 
pathways.show <- cellchat.ctrl@netP$pathways
groupSize <- as.numeric(table(cellchat.ctrl@idents)) # number of cells in each cell group
vertex.receiver = seq(1,9) # a numeric vector

#P1-- Hierarchy plot
netVisual_aggregate(cellchat.ctrl, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
#P2-- Circle plot
netVisual_aggregate(cellchat.ctrl, signaling = pathways.show, layout = "circle", vertex.size = groupSize)

#P3--Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
netAnalysis_contribution(cellchat.ctrl, signaling = pathways.show)

#Identify signaling roles of cell groups
cellchat.ctrl <- netAnalysis_signalingRole(cellchat.ctrl, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#P4--Visualize the signaling roles of cell groups
netVisual_signalingRole(cellchat.ctrl, signaling = pathways.show)

#Identify and visualize outgoing communication pattern of secreting cells
nPatterns = 8
cellchat.ctrl <- identifyCommunicationPatterns(cellchat.ctrl, pattern = "outgoing", k = nPatterns)
#P5-- river plot
netAnalysis_river(cellchat.ctrl, pattern = "outgoing")
#P6-- dot plot
netAnalysis_dot(cellchat.ctrl, pattern = "outgoing")
#Identify and visualize incoming communication pattern of target cells
nPatterns = 8
cellchat.ctrl <- identifyCommunicationPatterns(cellchat.ctrl, pattern = "incoming", k = nPatterns)
#P7 -- river plot
netAnalysis_river(cellchat.ctrl, pattern = "incoming")
#P8--- dot plot
netAnalysis_dot(cellchat.ctrl, pattern = "incoming")


##Manifold and classification learning analysis of signaling networks
#Identify signaling groups based on their functional similarity
cellchat.ctrl <- computeNetSimilarity(cellchat.ctrl, type = "functional")
cellchat.ctrl <- netEmbedding(cellchat.ctrl, type = "functional")
cellchat.ctrl <- netClustering(cellchat.ctrl, type = "functional")
#P9-- Visualization in 2D-space
netVisual_embedding(cellchat.ctrl, type = "functional")
netVisual_embeddingZoomIn(cellchat.ctrl, type = "functional")

#Identify signaling groups based on structure similarity
cellchat.ctrl <- computeNetSimilarity(cellchat.ctrl, type = "structural")
cellchat.ctrl <- netEmbedding(cellchat.ctrl, type = "structural")
cellchat.ctrl <- netClustering(cellchat.ctrl, type = "structural")
#P10-- Visualization in 2D-space
netVisual_embedding(cellchat.ctrl, type = "structural")
netVisual_embeddingZoomIn(cellchat.ctrl, type = "structural")

#P11-Bubble plot (to cells in x)
cellchat.ctrl<-rankNetPairwise(cellchat.ctrl)
pathways.show <- c('HGF') 
netVisual_bubble(cellchat.ctrl, from =c(4), to =(11))


#Select which database to use:
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") # use Secreted Signaling for cell-cell communication analysis
cellchat.lig@DB <- CellChatDB.use # set the used database in the object

#Preprocessing the expression data for cell-cell communication analysis
cellchat.lig <- subsetData(cellchat.lig) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat.lig <- identifyOverExpressedGenes(cellchat.lig)
cellchat.lig <- identifyOverExpressedInteractions(cellchat.lig)
cellchat.lig <- projectData(cellchat.lig, PPI.mouse)

#Inference of cell-cell communication network ##
#Compute the communication probability and infer cellular communication network
cellchat.lig <- computeCommunProb(cellchat.lig)
#fer the cell-cell communication at a signaling pathway level
cellchat.lig <- computeCommunProbPathway(cellchat.lig)
#Calculate the aggregated cell-cell communication network
cellchat.lig <- aggregateNet(cellchat.lig)

######Visualization and systems analysis of cell-cell communication network
#Visualize each signaling pathway using hierarchy plot or circle plot

pathways.show <- c('HGF') 
pathways.show <- cellchat.lig@netP$pathways
groupSize <- as.numeric(table(cellchat.lig@idents)) # number of cells in each cell group
vertex.receiver = seq(1,9) # a numeric vector


#P1-- Hierarchy plot
netVisual_aggregate(cellchat.lig, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
#P2-- Circle plot
netVisual_aggregate(cellchat.lig, signaling = pathways.show, layout = "circle", vertex.size = groupSize)

#P3--Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
netAnalysis_contribution(cellchat.lig, signaling = pathways.show)

#Identify signaling roles of cell groups
cellchat.lig <- netAnalysis_signalingRole(cellchat.lig, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#P4--Visualize the signaling roles of cell groups
netVisual_signalingRole(cellchat.lig, signaling = pathways.show)

#Identify and visualize outgoing communication pattern of secreting cells
nPatterns = 8
cellchat.lig <- identifyCommunicationPatterns(cellchat.lig, pattern = "outgoing", k = nPatterns)
#P5-- river plot
netAnalysis_river(cellchat.lig, pattern = "outgoing")
ggsave('.png',dpi=300,width = 4,height = 5,units = 'in')
#P6-- dot plot
netAnalysis_dot(cellchat.lig, pattern = "outgoing")
ggsave('.png',dpi=300,width = 8,height = 4,units = 'in')
#Identify and visualize incoming communication pattern of target cells
nPatterns = 8
cellchat.lig <- identifyCommunicationPatterns(cellchat.lig, pattern = "incoming", k = nPatterns)
#P7 -- river plot
netAnalysis_river(cellchat.lig, pattern = "incoming")
#P8--- dot plot
netAnalysis_dot(cellchat.lig, pattern = "incoming")


##Manifold and classification learning analysis of signaling networks
#Identify signaling groups based on their functional similarity
cellchat.lig <- computeNetSimilarity(cellchat.lig, type = "functional")
cellchat.lig <- netEmbedding(cellchat.lig, type = "functional")
cellchat.lig <- netClustering(cellchat.lig, type = "functional")
#P9-- Visualization in 2D-space
netVisual_embedding(cellchat.lig, type = "functional")
netVisual_embeddingZoomIn(cellchat.lig, type = "functional")

#Identify signaling groups based on structure similarity
cellchat.lig <- computeNetSimilarity(cellchat.lig, type = "structural")
cellchat.lig <- netEmbedding(cellchat.lig, type = "strucßtural")
cellchat.lig <- netClustering(cellchat.lig, type = "structural")
#P10-- Visualization in 2D-space
netVisual_embedding(cellchat.lig, type = "structural")
netVisual_embeddingZoomIn(cellchat.lig, type = "structural")

#P11-Bubble plot
cellchat.lig<-rankNetPairwise(cellchat.lig)
pathways.show <- c('VISTA') 
netVisual_bubble(cellchat.lig, from =c('Fibroblasts','Telocytes'), to =c('B cells','T cells'))
netVisual_bubble(cellchat.lig, from =c(1,2), to =c(5,6))


######Combine two cellchats####

cellchat.ctrl<-readRDS("cellchat_ctrl_contact.rds")
cellchat.lig<- readRDS('cellchat_lig_Cellcell_contact.rds')
cellchat <- mergeCellChat(list(cellchat.ctrl, cellchat.lig), add.names = c("Ctrl","lig"))

#Run manifold and classification learning analysis
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural") #4mins
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization
netVisual_embeddingPairwise(cellchat, type = "structural")
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural")

#Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, type = "structural")

#Identify and visualize the conserved and context-specific signaling pathways
rankNet(cellchat, mode = "comparison")


######Further Analysis of Cellchat######
#extract LR pairs
pairLR <- searchPair(signaling = "HGF", pairLR.use = CellChatDB$interaction, key = "pathway_name", matching.exact = T, pair.only = F)

####for comapring DEG of two datasets####
s.combined$compare <- paste(Idents(s.combined), s.combined$batch, sep = "_")
Idents(s.combined) <- "compare"
DefaultAssay(s.combined) <- "RNA"
m.response <- FindMarkers(s.combined, ident.1 = "Telocytes_Ligature", ident.2 = "Telocytes_Ctrl")
m.response$gene<-rownames(m.response)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
Idents(s.combined) <- 'my_clusters'
m <- subset(s.combined, idents = 'Telocytes')
Idents(m) <- "batch"
avg.m <- log1p(AverageExpression(m, verbose = FALSE)$RNA)
avg.m$gene <- rownames(avg.m)
genes.table = m.response %>%top_n(-10, avg_logFC)
genes.tablel = m.response %>%top_n(10, avg_logFC)
genes.to.label <-c(genes.table$gene,genes.tablel$gene)
p1<-ggplot(avg.m, aes(Ctrl, Ligature)) + geom_point() + ggtitle("Telocytes")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE,xnudge = 0,ynudge = 0)



######## convert to loom file ##############

library(loomR)
sdata<- subset(s.combined, idents = c('Macrophages 1','Macrophages 2','Macrophages 3','Macrophages 4')) #s is combined dataset
m<-FindAllMarkers(sdata)
marker <- m %>% group_by(cluster) %>% top_n(20, avg_log2FC)

sdata@graphs <- list() ####don't save Seurat after this, clean the enviroment and read RDS for further analysis!!
# seurat convert to loom file
sdata.loom <- as.loom(x = sdata, filename = "~/data/rmc1.loom", verbose = FALSE)
# Always remember to close loom files when done
sdata.loom$close_all()
# Jupyter: jupyter-notebook --no-browser --port=9999