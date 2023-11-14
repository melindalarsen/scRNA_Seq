
if (!require("pacman")) install.packages("pacman")
pacman::p_load(SingleCellExperiment, celda, Seurat, ggplot2, cowplot, scater, patchwork, 
               dplyr, sctransform, glmGamPoi, stringr, scCustomize, scDblFinder, singleCellTK, biomaRt,
               clustree, hdf5r, rhdf5, sctransform, remotes)

devtools::install_github("sqjin/CellChat")
library("CellChat")

dataPath <- "Z:/Kenney/E034 Amber reviewer comments/FinalDataSetsSCTransformIntegrated"
projectPath <- "Z:/Kenney/E034 Amber reviewer comments/JDRFigures"

setwd(projectPath)
set.seed(12354)
options(future.globals.maxSize = 10 * 1024^3)
indexList <- list(mock = 1, ligated = 2, deligated = 3, immune = 1, endothelial = 2, fibroblast = 3, 
                  epithelial = 4, pericyte = 5, glia = 6)

cellChatFromSeurat <- function(aSeurat, protocol){
  if(0)
  {
    aSeurat = aSeuratIntegrated
    protocol <- "All"
  }
  print(aSeurat)
  print(sprintf("Createing cellChat for %s", protocol))
  data.input <- GetAssayData(aSeurat, assay = "SCT", slot = "data")
  labels <- Idents(aSeurat)
  meta <- data.frame(group = labels, row.names = names(labels))
  
  acellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  
  acellchat <- addMeta(acellchat, meta = meta) 
  acellchat <- setIdent(acellchat, ident.use = "group") # set "group" as default cell identity
  levels(acellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(acellchat@idents)) # number of cells in each cell group
  
  CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
  #showDatabaseCategory(CellChatDB)
  
  # Show the structure of the database
  #dplyr::glimpse(CellChatDB$interaction)
  
  # use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  
  #Fixes random error that cellchat has with its own DB
  #https://github.com/sqjin/CellChat/issues/45
  removeID <- which(CellChatDB.use[["interaction"]]$ligand == "H2-BI") # 1887
  CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1887,]
  removeID <- which(CellChatDB.use[["interaction"]]$ligand == "H2-Ea-ps") #1900
  CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-removeID,]
  
  acellchat@DB <- CellChatDB.use
  
  # subset the expression data of signaling genes for saving computation cost
  acellchat <- subsetData(acellchat) # This step is necessary even if using the whole database
  future::plan("multisession", workers = 8) # do parallel
  
  acellchat <- identifyOverExpressedGenes(acellchat)
  acellchat <- identifyOverExpressedInteractions(acellchat)
  #acellchat <- projectData(acellchat, PPI.mouse)   #possable cleanup step, optional, need to use raw.use = FALSE in computeComunProb()
  
  acellchat <- computeCommunProb(acellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  acellchat <- filterCommunication(acellchat, min.cells = 10)
  
  acellchat <- computeCommunProbPathway(acellchat)
  
  acellchat <- aggregateNet(acellchat)
  
  # Compute the network centrality scores
  # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  acellchat <- netAnalysis_computeCentrality(acellchat, slot.name = "netP") 
  #Sets 'communication patterns'
  acellchat <- identifyCommunicationPatterns(acellchat, pattern = "outgoing", k = 3)
  acellchat <- identifyCommunicationPatterns(acellchat, pattern = "incoming", k = 3)
  
  acellchat@meta$protocol <- protocol
  
  return(acellchat)
}
saveClusterMarkers <- function(seuratObject, projName = "project"){
  
  #Debug assignment
  if(0){
    seuratObject = aSeurat
    projName = "cellbenderSeurat"
  }
  
  dir.create(projName)
  
  apath <- sprintf("%s/geneAverages_%s.csv", projName, projName)
  write.csv(FindAllMarkers(seuratObject, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
  apath <- sprintf("%s/positveGeneAverages_%s.csv", projName, projName)           
  write.csv(FindAllMarkers(seuratObject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
  apath <- sprintf("%s/cellsPerCluster%s.csv", projName, projName)
  write.csv(table(Idents(seuratObject)), file = apath)
  
  # colnames(tmp)[colnames(tmp) == "pct.1"] = "First"
  # colnames(tmp)[colnames(tmp) == "pct.2"] = "Second"
}
netVisual_bubble_MLD <- function(aCellChat, source = c(1), target = c(1), apath = "netVisualBubble") {
  
  netVisual_bubble_invert(aCellChat, source, target, comparison = c(indexList$mock, indexList$ligated),
                          label = c("Mock", "Ligated"), apath = apath)
  netVisual_bubble_invert(aCellChat, source, target, comparison = c(indexList$mock, indexList$deligated),
                          label = c("Mock", "Deligated"), apath = apath)
  netVisual_bubble_invert(aCellChat, source, target, comparison = c(indexList$deligated, indexList$ligated),
                          label = c("Deligated", "Ligated"), apath = apath)
}
netVisual_bubble_invert <- function(aCellChat, source = c(1), target = c(1), comparison, label, apath = "netVisualBubble") {
  dir.create(apath)
  dir.create(sprintf("%s/To_%s", apath, levels(aCellChat@meta$group)[source[1]]))
  dir.create(sprintf("%s/To_%s", apath, levels(aCellChat@meta$group)[target[1]]))
  
  fileNamePath <- sprintf("%s/To_%s/netVisual_bubble_%sv%s_From_%s_to_%s.pdf", apath, 
                          levels(aCellChat@meta$group)[target[1]], label[1], label[2], 
                          levels(aCellChat@meta$group)[source[1]], levels(aCellChat@meta$group)[target[1]])
  pdf(fileNamePath)
  tryCatch({
    print(netVisual_bubble(aCellChat, sources.use = source, targets.use = target,  
                           comparison = comparison, 
                           max.dataset = comparison[1], 
                           title.name = sprintf("Increased signaling in %s vs %s", label[1], label[2]), 
                           angle.x = 45, 
                           remove.isolate = T))
  }, error=function(e){cat(fileNamePath, conditionMessage(e), "\n")})
  dev.off()
  
  fileNamePath <- sprintf("%s/To_%s/netVisual_bubble_%sv%s_From_%s_to_%s.pdf", apath, 
                          levels(aCellChat@meta$group)[target[1]], label[2], label[1], 
                          levels(aCellChat@meta$group)[source[1]], levels(aCellChat@meta$group)[target[1]])
  pdf(fileNamePath)
  tryCatch({
    print(netVisual_bubble(aCellChat, sources.use = source, targets.use = target,  
                           comparison = comparison, 
                           max.dataset = comparison[2], 
                           title.name = sprintf("Increased signaling in %s vs %s", label[2], label[1]), 
                           angle.x = 45, 
                           remove.isolate = T))
  }, error=function(e){cat(fileNamePath ,conditionMessage(e), "\n")})
  dev.off()
}
convert_human_to_mouse <- function(gene_list) {
  #Adapted from https://www.biostars.org/p/9567892/
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (output[,2])
}

#Load seurat object of Integrated datasets created in E034 and clutered with a resolution of 2
aSeuratIntegrated <- readRDS(sprintf("%s/AmberAllSCTransformIntegrated-2.rds", dataPath))

Idents(aSeuratIntegrated) <- "seurat_clusters"
levels(aSeuratIntegrated)
pdf("UMAP_integrate_All.pdf", width = 40, height = 30)
print(DimPlot(aSeuratIntegrated, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()

Idents(aSeuratIntegrated) <- "protocol"
levels(aSeuratIntegrated)
aSeuratIntegrated$protocol <- factor(aSeuratIntegrated$protocol, levels = c("Mock", "Ligated", "Deligated"))
Idents(aSeuratIntegrated) <- "protocol"
levels(aSeuratIntegrated)
pdf("UMAP_integrate_All_protocol.pdf", width = 40, height = 30)
print(DimPlot(aSeuratIntegrated, reduction = "umap", label = FALSE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 120)))
dev.off()

#Rename cluster IDs to cell type
Idents(aSeuratIntegrated) <- "seurat_clusters"
levels(aSeuratIntegrated)
newClusterIDsLong <- c("0_Immune", "1_Endothelial", "2_Endothelial", "3_Endothelial", "4_Endothelial", "5_Immune", 
                       "6_Endothelial", "7_Immune", "8_Immune", "9_Fibroblast", "10_Endothelial", 
                       "11_Endothelial", "12_Immune", "13_Endothelial", "14_Immune", "15_Endothelial", 
                       "16_Endothelial", "17_Endothelial", "18_Epithelial", "19_Fibroblast", "20_Pericyte", 
                       "21_Endothelial", "22_Immune", "23_Endothelial", "24_Pericyte", "25_Pericyte",
                       "26_Epithelial", "27_Immune", "28_Endothelial", "29_Glia", "30_Immune",
                       "31_Epithelial", "32_Immune", "33_Immune", "34_Endothelial", "35_Endothelial",
                       "36_Endothelial", "37_Immune", "38_Endothelial", "39_Epithelial", "40_Endothelial",
                       "41_Immune", "42_Immune", "43_Immune", "44_Fibroblast", "45_Immune",
                       "46_Fibroblast", "47_Endothelial", "48_Immune")
newClusterIDsShort <- c("Immune", "Endothelial", "Endothelial", "Endothelial", "Endothelial", "Immune", 
                        "Endothelial", "Immune", "Immune", "Fibroblast", "Endothelial", 
                        "Endothelial", "Immune", "Endothelial", "Immune", "Endothelial", 
                        "Endothelial", "Endothelial", "Epithelial", "Fibroblast", "Pericyte", 
                        "Endothelial", "Immune", "Endothelial", "Pericyte", "Pericyte",
                        "Epithelial", "Immune", "Endothelial", "Glia", "Immune",
                        "Epithelial", "Immune", "Immune", "Endothelial", "Endothelial",
                        "Endothelial", "Immune", "Endothelial", "Epithelial", "Endothelial",
                        "Immune", "Immune", "Immune", "Fibroblast", "Immune",
                        "Fibroblast", "Endothelial", "Immune")
Idents(aSeuratIntegrated) <- "seurat_clusters"
newClusterIDs <- newClusterIDsShort
names(newClusterIDs) <- levels(aSeuratIntegrated)
aSeuratIntegrated <- RenameIdents(aSeuratIntegrated, newClusterIDs)
levels(aSeuratIntegrated)

pdf("UMAP_integrate_All_NewID.pdf", width = 40, height = 30)
print(DimPlot(aSeuratIntegrated, reduction = "umap", label = TRUE, label.size = 40, pt.size = 0.5, repel = TRUE) +
        NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()

#Generate the files for NATMI, NATMI needs to be run on the commandline
seuratList <- c()
seuratList <- c(seuratList, subset(aSeuratIntegrated, protocol == "Mock"))
seuratList <- c(seuratList, subset(aSeuratIntegrated, protocol == "Ligated"))
seuratList <- c(seuratList, subset(aSeuratIntegrated, protocol == "Deligated"))

table(Idents(seuratList[[indexList$mock]]))[2]
ncol(seuratList[[indexList$mock]])
table(Idents(seuratList[[indexList$ligated]]))[2]
ncol(seuratList[[indexList$ligated]])
table(Idents(seuratList[[indexList$deligated]]))[2]
ncol(seuratList[[indexList$deligated]])
      
for(i in seuratList){
  #From NATMI
  write.csv(100 * (exp(as.matrix(GetAssayData(object = i, assay = "SCT", slot = "data"))) - 1), 
            sprintf("data-%s.csv", i@meta.data$protocol[1]), row.names = T)
  write.csv(Idents(object = i), sprintf("metadata-%s.csv", i@meta.data$protocol[1]), row.names = T)
  
}

cellchatList <- c()
for(i in 1:3){
  aCellChat <- cellChatFromSeurat(seuratList[[i]], seuratList[[i]]$protocol[1])
  saveRDS(aCellChat, file = sprintf("CellChat%s.rds", seuratList[[i]]$protocol[1]))
  cellchatList <- c(cellchatList, aCellChat)
}

aCellchatMerge <- mergeCellChat(cellchatList, cell.prefix = TRUE, add.names = c(cellchatList[[1]]@meta$protocol[1], 
                                                            cellchatList[[2]]@meta$protocol[1], 
                                                            cellchatList[[3]]@meta$protocol[1]))
saveRDS(aCellchatMerge, file = "CellChatMerge.rds")

aCellchatMerge <- computeNetSimilarityPairwise(aCellchatMerge, type = "structural")
aCellchatMerge <- netEmbedding(aCellchatMerge, type = "structural", umap.method = "uwot")
aCellchatMerge <- netClustering(aCellchatMerge, type = "structural", do.parallel = FALSE)

aCellchatMerge <- computeNetSimilarityPairwise(aCellchatMerge, type = "functional")
aCellchatMerge <- netEmbedding(aCellchatMerge, type = "functional", umap.method = "uwot")
aCellchatMerge <- netClustering(aCellchatMerge, type = "functional", do.parallel = FALSE)

saveRDS(aCellchatMerge, file = "CellChatMerge2.rds")

cellchatList <- c()
aCellChat <- readRDS("CellChatMock.rds")
cellchatList <- c(cellchatList, aCellChat)
aCellChat <- readRDS("CellChatLigated.rds")
cellchatList <- c(cellchatList, aCellChat)
aCellChat <- readRDS("CellChatDeligated.rds")
cellchatList <- c(cellchatList, aCellChat)
aCellchatMerge <- readRDS("CellChatMerge2.rds")

levels(cellchatList[[1]]@idents)

#Differanl gene expression---------------------
levels(aSeuratIntegrated)

aSeuratIntegrated$protocol_cluster <- paste0(aSeuratIntegrated$protocol, "_", aSeuratIntegrated@active.ident)
Idents(aSeuratIntegrated) <- "protocol_cluster"
levels(aSeuratIntegrated)

apath <- "Pos_Markers_Mock_Endo_Ligated_Endo.csv"
write.csv(FindMarkers(aSeuratIntegrated, ident.1 = "Mock_Endothelial", ident.2 = "Ligated_Endothelial", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
apath <- "Pos_Markers_Ligated_Endo_Mock_Endo.csv"
write.csv(FindMarkers(aSeuratIntegrated, ident.1 = "Ligated_Endothelial", ident.2 = "Mock_Endothelial", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
apath <- "Pos_Markers_Mock_Endo_Deligated_Endo.csv"
write.csv(FindMarkers(aSeuratIntegrated, ident.1 = "Mock_Endothelial", ident.2 = "Deligated_Endothelial", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
apath <- "Pos_Markers_Deligated_Endo_Mock_Endo.csv"
write.csv(FindMarkers(aSeuratIntegrated, ident.1 = "Deligated_Endothelial", ident.2 = "Mock_Endothelial", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
apath <- "Pos_Markers_Ligated_Endo_Deligated_Endo.csv"
write.csv(FindMarkers(aSeuratIntegrated, ident.1 = "Ligated_Endothelial", ident.2 = "Deligated_Endothelial", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
apath <- "Pos_Markers_Deligated_Endo_Ligated_Endo.csv"
write.csv(FindMarkers(aSeuratIntegrated, ident.1 = "Deligated_Endothelial", ident.2 = "Ligated_Endothelial", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)


#Netvisual bubble show signaling increase from cluster to cluster between conditions (Mock, Ligated, Deligated)---------------------------------
#From Endotheal cells
levels(aCellchatMerge@meta$group)


for(i in 1:nlevels(aCellchatMerge@meta$group)){
  netVisual_bubble_MLD(aCellchatMerge, source = 2, target = i, apath = "netVisualBubble")
}

#General cellchat figures---------------------------

groupSize <- as.numeric(table(aCellChat@idents))

pdf("netVisual_circle_End_weight_mock.pdf")
netVisual_circle(cellchatList[[indexList$mock]]@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use = "Endothelial")
dev.off()

pdf("netVisual_circle_End_weight_ligated.pdf")
netVisual_circle(cellchatList[[indexList$ligated]]@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use = "Endothelial")
dev.off()

pdf("netVisual_circle_End_weight_deligated.pdf")
netVisual_circle(cellchatList[[indexList$deligated]]@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use = "Endothelial")
dev.off()

pdf("netVisual_circle_End_count_mock.pdf")
netVisual_circle(cellchatList[[indexList$mock]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use = "Endothelial")
dev.off()

pdf("netVisual_circle_End_count_ligated.pdf")
netVisual_circle(cellchatList[[indexList$ligated]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use = "Endothelial")
dev.off()

pdf("netVisual_circle_End_count_deligated.pdf")
netVisual_circle(cellchatList[[indexList$deligated]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use = "Endothelial")
dev.off()

pathway.union <- union(cellchatList[[indexList$mock]]@netP$pathways, 
                       cellchatList[[indexList$ligated]]@netP$pathways)
pathway.union <- union(pathway.union, 
                       cellchatList[[indexList$deligated]]@netP$pathways)

pdf("netAnalysis_signalingRole_heatmap_mock.pdf", height = 9)
  print(netAnalysis_signalingRole_heatmap(cellchatList[[indexList$mock]], pattern = "all", 
                                        signaling = pathway.union, title = "Mock", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
pdf("netAnalysis_signalingRole_heatmap_ligated.pdf", height = 9)
  print(netAnalysis_signalingRole_heatmap(cellchatList[[indexList$ligated]], pattern = "all", 
                                        signaling = pathway.union, title = "Ligated", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
pdf("netAnalysis_signalingRole_heatmap_deligated.pdf", height = 9)
  print(netAnalysis_signalingRole_heatmap(cellchatList[[indexList$deligated]], pattern = "all", 
                                        signaling = pathway.union, title = "Deligated", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()

pdf("compareInteractions_count.pdf")
compareInteractions(aCellchatMerge, show.legend = F, group = c(1,2,3), width = .2)
dev.off()

pdf("compareInteractions_strength.pdf")
compareInteractions(aCellchatMerge, show.legend = F, group = c(1,2,3), measure = "weight", width = .2, digits=2)
dev.off()

#reordering the clusters so that netVisual_diffInteraction works. netVisual_diffInteraction is bugged and will only work for clusters 1 and <max>
tmp <- c()
aCellChat <- updateClusterLabels(cellchatList[[1]], new.order = c("Endothelial", "Immune", "Fibroblast", "Epithelial", "Pericyte", "Glia"))
tmp <- c(tmp, aCellChat)
aCellChat <- updateClusterLabels(cellchatList[[2]], new.order = c("Endothelial", "Immune", "Fibroblast", "Epithelial", "Pericyte", "Glia"))
tmp <- c(tmp, aCellChat)
aCellChat <- updateClusterLabels(cellchatList[[3]], new.order = c("Endothelial", "Immune", "Fibroblast", "Epithelial", "Pericyte", "Glia"))
tmp <- c(tmp, aCellChat)

Mergetmp <- mergeCellChat(tmp, cell.prefix = TRUE, add.names = c(tmp[[1]]@meta$protocol[1], 
                                                                                tmp[[2]]@meta$protocol[1], 
                                                                                tmp[[3]]@meta$protocol[1]))

Mergetmp <- computeNetSimilarityPairwise(Mergetmp, type = "structural")
Mergetmp <- netEmbedding(Mergetmp, type = "structural", umap.method = "uwot")
Mergetmp <- netClustering(Mergetmp, type = "structural", do.parallel = FALSE)

Mergetmp <- computeNetSimilarityPairwise(Mergetmp, type = "functional")
Mergetmp <- netEmbedding(Mergetmp, type = "functional", umap.method = "uwot")
Mergetmp <- netClustering(Mergetmp, type = "functional", do.parallel = FALSE)

pdf("cellchat/netVisual_diffInteraction_mockvligated.pdf")
netVisual_diffInteraction(Mergetmp, comparison = c(indexList$mock, indexList$ligated), measure = "weight",
                          sources.use = 1)
dev.off()
pdf("cellchat/netVisual_diffInteraction_ligatedvmock.pdf")
netVisual_diffInteraction(Mergetmp, comparison = c(indexList$ligated, indexList$mock), measure = "weight",
                          sources.use = 1)
dev.off()
pdf("cellchat/netVisual_diffInteraction_mockvdeligated.pdf")
netVisual_diffInteraction(Mergetmp, comparison = c(indexList$mock, indexList$deligated), measure = "weight",
                          sources.use = 1)
dev.off()
pdf("cellchat/netVisual_diffInteraction_deligatedvmock.pdf")
netVisual_diffInteraction(Mergetmp, comparison = c(indexList$deligated, indexList$mock), measure = "weight",
                          sources.use = 1)
dev.off()
pdf("cellchat/netVisual_diffInteraction_ligatedvdeligated.pdf")
netVisual_diffInteraction(Mergetmp, comparison = c(indexList$ligated, indexList$deligated), measure = "weight",
                          sources.use = 1)
dev.off()
pdf("cellchat/netVisual_diffInteraction_deligatedvligated.pdf")
netVisual_diffInteraction(Mergetmp, comparison = c(indexList$deligated, indexList$ligated), measure = "weight",
                          sources.use = 1)
dev.off()

#EndMT and Nolan factors dotplots------------------------------------

netVisual_diffInteraction(aCellchatMerge, comparison = c(indexList$mock, indexList$ligated),
                           sources.use = c("Endothelial"))

NolanList <- c("Bmp2", "Ccl2", "Ccl3", "Ccl6", "Ccl9", "Col14a1", "Col4a1", "Cxcl12", "Cxcl13", "Cxcl2", "Dkk2", 
               "Edn1", "Egfl7", "Eln", "Esm1", "Fgf7", "Fn1", "Gja1", "Gja4", "Igf1", "Il1a", "Il33", "Jag1", 
               "Kitl", "Lama4", "Lamb1", "Lamc1", "Mmp13", "Mmp9", "Pdgfa", "Tgfb2", "Tgfb3", "Tnf", "Vegfc")

pdf("Nolan-Endo-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = NolanList, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 20, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 70), axis.text = element_text(size = 70), 
              legend.text = element_text(size = 70), legend.key.size = unit(2, "cm"), 
              legend.title  = element_text(size = 70), legend.spacing = unit(2, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

#ENdoMT list(s)
#Emerging roles of inflammation-mediated endothelial–mesenchymal transition in health and disease
#Yasuhiro Yoshimatsu & Tetsuro Watabe (2022)
EndoMT_Inducer <- c("Bmp9", "Mmp14", "Tgfb1", "Tgfb2")
EndoMT_Inflammation <- c("Ccl2", "Ccl3", "Cxcl2", "Icam1", "Ifng", "Il1b", "Il4", "Il6", "Jak2", "Rela", "Stat3", 
                         "Tnf", "Vcam1", "C1qa", "C1qb", "C5ar1")
EndoMT_EC_Gene <- c("Apln", "Cdh5", "Erg", "Fgfr1", "Flt4", "Fstl3", "Kdr", "H19", "Lyve1", 
                    "Nos3", "Pdpn", "Pecam1", "Prox1", "Tal1", "Tie1", 
                    "Tek", "Vwf", "Wnt5a")

EndoMT_Mesen_Gene <- c("Acta2", "Cdh2", "Cnn1", "Col12a1", "Col1a1", "Col1a2", "Col3a1", 
                       "Col6a1", "Ctgf", "Fap", "Fn1", "Lamc2", "Ly6a", "Atxn1", "Mmp14", "Mmp2", "Mmp9", "Notch3", 
                       "Pcolce", "Pdgfrb", "S100a4", "Aifm2", "Serpine1", "Tagln", "Tgfbi", "Tnc", "Vim")

EndoMT_Mesen_Gene_short <- c("Acta2", "Col12a1", "Col1a1", "Col1a2", "Col3a1", 
                       "Col6a1", "Fn1", "Lamc2", "Ly6a", "Atxn1", "Notch3", 
                       "Pdgfrb", "S100a4", "Vim")

EndoMT_TF <- c("Hes1", "Hey1", "Hey2", "Mkl1", "Mrtfa", "Snai1", "Snai2", "Twist1", "Twist2", "Zeb1", "Zeb2")

EndoMT_subset <- c("Tgfb1", "Tgfb2",  #Inducer
                   "Icam1", "Jak2", "Rela", "Stat3", #Inflammation
                   "Cdh5", "Erg", "Flt4", "Kdr", "Nos3", "Pecam1", "Tal1", "Tie1", "Tek", #EC Genes
                   "Ly6a", "Atxn1", "Vim", #Mes Genes
                   "Hes1", "Hey1", "Snai1", "Snai2", "Zeb1") #Transcription factors

pdf("EndoMT-Endo-Inducer-DotPlot.pdf", width = 40, height = 30)
print(DotPlot(aSeuratIntegrated, features = EndoMT_Inducer, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 60, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 60), axis.text = element_text(size = 50), 
              legend.text = element_text(size = 40), legend.key.size = unit(3, 'cm'), 
              legend.title  = element_text(size = 40), legend.spacing = unit(3, 'cm')) +
        RotatedAxis())
dev.off()

pdf("EndoMT-Endo-Inflammation-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = EndoMT_Inflammation, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 60, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"), 
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

pdf("EndoMT-Endo-EC_Gene-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = EndoMT_EC_Gene, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 60, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 80), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 60), legend.key.size = unit(2, "cm"), 
              legend.title  = element_text(size = 60), legend.spacing = unit(2, "cm"), 
              axis.text.x = element_text(angle=90)))
dev.off()

pdf("EndoMT-Endo-Mesen_Gene_short-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = EndoMT_Mesen_Gene_short, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 30, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"), 
              axis.text.x = element_text(angle=90)))
dev.off()

pdf("EndoMT-Endo-Mesen_Gene-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = EndoMT_Mesen_Gene, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 30, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 50), axis.text = element_text(size = 90), 
              legend.text = element_text(size = 40), legend.key.size = unit(1, "cm"), 
              legend.title  = element_text(size = 50), legend.spacing = unit(1, "cm"), 
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

pdf("EndoMT-Endo-TF_Gene-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = EndoMT_TF, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 60, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

pdf("EndoMT-Endo-Subset-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = EndoMT_EC_Gene, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 40, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

#Endothelial subsetting---------------------------

aSeuratEndo <- subset(aSeuratIntegrated, idents = "Endothelial")
DefaultAssay(aSeuratEndo) <- "integrated"

aSeuratEndo <- SCTransform(aSeuratEndo, vst.flavor = "v2", verbose = FALSE)
aSeuratEndo <- RunPCA(aSeuratEndo, npcs = 40, verbose = FALSE)
aSeuratEndo <- RunUMAP(aSeuratEndo, reduction = "pca", dims = 1:40)
aSeuratEndo <- FindNeighbors(aSeuratEndo, reduction = "pca", dims = 1:40)

clusterResolution <- 0.5
seuratTree <- FindClusters(aSeuratEndo, resolution = (seq(0, (clusterResolution*2), by = (clusterResolution/5))))

pdf("Endo-clusttree.pdf", width = 40, height = 30)
print(clustree(seuratTree, prefix = "SCT_snn_res."))
dev.off()

aSeuratEndo <- FindClusters(aSeuratEndo, resolution = clusterResolution)
aSeuratEndo <- RunTSNE(aSeuratEndo, dims = 1:40, check_duplicates=FALSE)
aSeuratEndo <- RunUMAP(aSeuratEndo, dims = 1:40)

Idents(aSeuratEndo) <- "seurat_clusters"
levels(aSeuratEndo)
pdf("Endo-umap.pdf", width = 40, height = 30)
print(DimPlot(aSeuratEndo, reduction = "umap", label = TRUE, label.size = 40, pt.size = 0.5, repel = TRUE) +
        NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40)))
dev.off()

DefaultAssay(aSeuratEndo) <- "SCT"
saveClusterMarkers(aSeuratEndo, projName = "Endo-markers")
  
Artery_Markers_Kaluka <- c("8430408G22Rik", "Clu", "Crip1", "Fbln2", "Gja4", "Hey1", "Mecom", 
                           "Sat1", "Sema3g", "Sox17", "Tm4sf1", "Tsc22d1")

Vein_Markers_Kaluka <- c("Apoe", "Bgn", "Ctla2a", "Icam1", "Il6st", "Ptgs1", "Tmsb10", "Vcam1", "Vwf")

Capillary_Markers_Kaluka <- c("AW112010", "BC028528", "Car4", "Cd200", "Cd300lg", "Gpihbp1", "Kdr", "Rgcc", 
                              "Sgk1", "Sparc")

Lymphatic_Markers_Kaluka <- c("Ccl21a", "Cd63", "Cp", "Fgl2", "Flt4", "Fth1", "Fxyd6", "Maf", "Marcks", "Mmrn1", 
                              "Pard6g", "Pdpn", "Prelp", "Reln", "Rnase4", "Scn1b", "Stab1", "Thy1", "Timp2", "Timp3")         

pdf("Markers_Endo_Artery.pdf", width = 40, height = 30)
print(DotPlot(aSeuratEndo, features = Artery_Markers_Kaluka, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 30) +
          theme(axis.title = element_text(size = 60), axis.text = element_text(size = 50), 
          legend.text = element_text(size = 40), legend.key.size = unit(3, 'cm'), 
          legend.title  = element_text(size = 40), legend.spacing = unit(3, 'cm')) +
          RotatedAxis())
dev.off()

pdf("Markers_Endo_Vein.pdf", width = 40, height = 30)
print(DotPlot(aSeuratEndo, features = Vein_Markers_Kaluka, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 30) +
        theme(axis.title = element_text(size = 60), axis.text = element_text(size = 50), 
              legend.text = element_text(size = 40), legend.key.size = unit(3, 'cm'), 
              legend.title  = element_text(size = 40), legend.spacing = unit(3, 'cm')) +
        RotatedAxis())
dev.off()

pdf("Markers_Endo_Capillary.pdf", width = 40, height = 30)
print(DotPlot(aSeuratEndo, features = Capillary_Markers_Kaluka, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 30) +
        theme(axis.title = element_text(size = 60), axis.text = element_text(size = 50), 
              legend.text = element_text(size = 40), legend.key.size = unit(3, 'cm'), 
              legend.title  = element_text(size = 40), legend.spacing = unit(3, 'cm')) +
        RotatedAxis())
dev.off()

pdf("Markers_Endo_Lymphatic.pdf", width = 40, height = 30)
print(DotPlot(aSeuratEndo, features = Lymphatic_Markers_Kaluka, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 30) +
        theme(axis.title = element_text(size = 60), axis.text = element_text(size = 50), 
              legend.text = element_text(size = 40), legend.key.size = unit(3, 'cm'), 
              legend.title  = element_text(size = 40), legend.spacing = unit(3, 'cm')) +
        RotatedAxis())
dev.off()

endo_lables_long_dot <- c("0_Capillary", "1_Capillary", "2_Capillary", "3_Capillary", "4_Artery", "5_Capillary", 
                      "6_Vein", "7_Artery", "8_Capillary", "9_Capillary", "10_Capillary", "11_Capillary", "12_RBC", 
                      "13_Capillary", "14_Capillary", "15_Lymphatic")

endo_lables_short_dot <- c("Capillary", "Capillary", "Capillary", "Capillary", "Artery", "Capillary", 
                          "Vein", "Artery", "Capillary", "Capillary", "Capillary", "Capillary", "RBC", 
                          "Capillary", "Capillary", "Lymphatic")

endo_lables_long <- c("0_Capillary", "1_Capillary", "2_Capillary", "3_Lymphatic", "4_Artery", "5_Lymphatic", 
                      "6_Vein", "7_Artery", "8_Endothelium", "9_Capillary", "10_Vein", "11_Capillary", "12_RBC", 
                      "13_Capillary", "14_Lymphatic", "15_Lymphatic")

endo_lables_short <- c("Capillary", "Capillary", "Capillary", "Lymphatic", "Artery", "Lymphatic", 
                      "Vein", "Artery", "Endothelium", "Capillary", "Vein", "Capillary", "RBC", 
                      "Capillary", "Lymphatic", "Lymphatic")



subsetIDMarkers<- c("Apoe", "Il6st", "Sema3g", "Clu", "Ccl21a", "Mmrn1", "Gpihbp1", "Rgcc", "Kdr", "Hbb-bs", "Hba-a1")
pdf("ClusterMarkers_subset_clusterID.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = subsetIDMarkers, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"), 
              axis.text.x = element_text(angle= 90, vjust = 0.45, hjust = 0.2)) +
        labs(x = "Vein Artery Lymph Capillary RBC"))
dev.off()

Idents(aSeuratEndo) <- "seurat_clusters"
newClusterIDs <- endo_lables_short_dot
names(newClusterIDs) <- levels(aSeuratEndo)
aSeuratEndo <- RenameIdents(aSeuratEndo, newClusterIDs)
levels(aSeuratEndo)
pdf("ClusterMarkers_subsetID.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = subsetIDMarkers, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"), 
              axis.text.x = element_text(angle= 90, vjust = 0.45)))
dev.off()

pdf("UMAP_Endo_NewID_dot.pdf", width = 40, height = 30)
print(DimPlot(aSeuratEndo, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()

#dotplot of markers after cluster lableing
aSeuratEndo$ID <- aSeuratEndo@active.ident
aSeuratEndo$protocol_subcluster <- paste0(aSeuratEndo$protocol, '_', aSeuratEndo$ID)
Idents(aSeuratEndo) <- "protocol_subcluster"
levels(aSeuratEndo)
table(Idents(aSeuratEndo))
levels(aSeuratEndo) <- c("Deligated_Vein", "Ligated_Vein", "Mock_Vein", 
                        "Deligated_Lymphatic", "Ligated_Lymphatic", "Mock_Lymphatic", 
                        "Deligated_Capillary", "Ligated_Capillary", "Mock_Capillary", 
                        "Deligated_Artery", "Ligated_Artery", "Mock_Artery", 
                        "Mock_RBC", "Ligated_RBC", "Deligated_RBC")

levels(aSeuratEndo)

endoSubsetIdents <- c("Deligated_Vein", "Ligated_Vein", "Mock_Vein", 
                      "Deligated_Lymphatic", "Ligated_Lymphatic", "Mock_Lymphatic", 
                      "Deligated_Capillary", "Ligated_Capillary", "Mock_Capillary", 
                      "Deligated_Artery", "Ligated_Artery", "Mock_Artery")

pdf("EmailedMarkers_subset.pdf", width = 40, height = 30)
print(DotPlot(aSeuratEndo, features = subsetIDMarkers, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30) +
                theme(axis.title = element_text(size = 60), axis.text = element_text(size = 50), 
                      legend.text = element_text(size = 40), legend.key.size = unit(3, 'cm'), 
                      legend.title  = element_text(size = 40), legend.spacing = unit(3, 'cm')) +
                RotatedAxis())
dev.off()

NolanListSubset <- c("Kitl", "Col14a1", "Egfl7", "Gja1", "Gja4", "Lama4", "Cxcl12")

pdf("./Endo-markers/NolanMarkers_subset_short.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = NolanListSubset, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30,
              idents = endoSubsetIdents) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

TFlist_short <- c("Hes1", "Snai1", "Snai2", "Zeb1")
pdf("./Endo-markers/TFMarkers_subset_short.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = TFlist_short, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30,
              idents = endoSubsetIdents) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

Meselist_short <- c("Ly6a", "Vim")
pdf("./Endo-markers/MeseMarkers_subset_short.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = Meselist_short, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30,
              idents = endoSubsetIdents) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

Endolist_short <- c("Cdh5", "Fgfr1", "Flt4", "Kdr", "Nos3", "Pecam1", "Tal1", "Tie1", "Tek", "Vwf")
pdf("./Endo-markers/EndoMarkers_subset_short.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = Endolist_short, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30,
              idents = endoSubsetIdents) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 80), legend.key.size = unit(2, "cm"), 
              legend.title  = element_text(size = 80), legend.spacing = unit(2, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

Infllist_short <- c("Jak2", "Rela", "Stat3")
pdf("./Endo-markers/InflMarkers_subset_short.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = Infllist_short, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30,
              idents = endoSubsetIdents) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()



table(Idents(aSeuratEndo))
table(Idents(subset(aSeuratEndo, protocol == "Mock")))
table(Idents(subset(aSeuratEndo, protocol == "Ligated")))
table(Idents(subset(aSeuratEndo, protocol == "Deligated")))

saveRDS(aSeuratEndo, file = "aSeuratEndo.rds")
aSeuratEndo <- readRDS("aSeuratEndo.rds")

DefaultAssay(aSeuratEndo)
levels(aSeuratEndo)
levels(aSeuratIntegrated)
aSeuratIntegrated$sub_cluster <- aSeuratIntegrated@active.ident

endoList <- c()
endoList <- c(endoList, subset(aSeuratEndo, idents = "Capillary"))
endoList <- c(endoList, subset(aSeuratEndo, idents = "Artery"))
endoList <- c(endoList, subset(aSeuratEndo, idents = "Vein"))
endoList <- c(endoList, subset(aSeuratEndo, idents = "RBC"))
endoList <- c(endoList, subset(aSeuratEndo, idents = "Lymphatic"))

aSeuratIntegrated$ID <- aSeuratIntegrated@active.ident
aSeuratEndo$ID <- aSeuratEndo@active.ident

levels(aSeuratIntegrated$ID)
levels(aSeuratEndo$ID)

aSeuratIntegratedNoEndo <- subset(aSeuratIntegrated, 
                                  ID == "Immune" | 
                                  ID == "Fibroblast" | 
                                  ID == "Epithelial" | 
                                  ID == "Pericyte" | 
                                  ID == "Glia")

aSeurat_subset <- merge(aSeuratIntegratedNoEndo, y = c(aSeuratEndo))

DefaultAssay(aSeurat_subset) <- "RNA"
aSeurat_subset <- SCTransform(aSeurat_subset, vst.flavor = "v2", verbose = FALSE)
DefaultAssay(aSeurat_subset) <- "SCT"

saveRDS(aSeurat_subset, file = "aSeuratSubset.rds")
aSeurat_subset <- readRDS("aSeuratSubset.rds")
levels(aSeurat_subset)

Cell_Cycle_Genes_of_Interest <- c("Cenpa", "Cenpe", "Cenpf", "Mki67",
                                  "Hist1h1a", "Hist1h1b", "Top2a", "Hist1h2ae", "Pcna", "Mcm2")
Idents(aSeuratEndo) <- "protocol_subcluster"
levels(aSeuratEndo)
pdf("./Endo-markers/Cell_Cycle_Endo_Subset.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = Cell_Cycle_Genes_of_Interest, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30,
              idents = endoSubsetIdents) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"), 
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

s.genes <- convert_human_to_mouse(cc.genes$g2m.genes)
g2m.genes <- convert_human_to_mouse(cc.genes$s.genes)

s.genes
g2m.genes

aSeuratEndo <- CellCycleScoring(aSeuratEndo, s.features = s.genes, g2m.features = g2m.genes)
levels(aSeuratEndo)
for(i in levels(aSeuratEndo)){
  tmp <- subset(aSeuratEndo, idents = i)
  Idents(tmp) <- "Phase"
  print(i)
  print(table(Idents(tmp)))
}

#Cellchat of sub-clustering--------------

aSeurat_subset <- readRDS("aSeuratSubset.rds")

seuratList <- c()
seuratList <- c(seuratList, subset(aSeurat_subset, protocol == "Mock"))
seuratList <- c(seuratList, subset(aSeurat_subset, protocol == "Ligated"))
seuratList <- c(seuratList, subset(aSeurat_subset, protocol == "Deligated"))

cellchatList <- c()
for(i in 1:3){
  aCellChat <- cellChatFromSeurat(seuratList[[i]], seuratList[[i]]$protocol[1])
  saveRDS(aCellChat, file = sprintf("CellChat_subset%s.rds", seuratList[[i]]$protocol[1]))
  cellchatList <- c(cellchatList, aCellChat)
}

aCellchatMerge <- mergeCellChat(cellchatList, cell.prefix = TRUE, add.names = c(cellchatList[[1]]@meta$protocol[1], 
                                                                                cellchatList[[2]]@meta$protocol[1], 
                                                                                cellchatList[[3]]@meta$protocol[1]))
saveRDS(aCellchatMerge, file = "CellChat_subsetMerge.rds")

aCellchatMerge <- computeNetSimilarityPairwise(aCellchatMerge, type = "structural")
aCellchatMerge <- netEmbedding(aCellchatMerge, type = "structural", umap.method = "uwot")
aCellchatMerge <- netClustering(aCellchatMerge, type = "structural", do.parallel = FALSE)

aCellchatMerge <- computeNetSimilarityPairwise(aCellchatMerge, type = "functional")
aCellchatMerge <- netEmbedding(aCellchatMerge, type = "functional", umap.method = "uwot")
aCellchatMerge <- netClustering(aCellchatMerge, type = "functional", do.parallel = FALSE)

saveRDS(aCellchatMerge, file = "CellChat_subsetMerge.rds")

cellchatSubList <- c()
aCellChat <- readRDS("CellChat_subsetMock.rds")
cellchatSubList <- c(cellchatSubList, aCellChat)
aCellChat <- readRDS("CellChat_subsetLigated.rds")
cellchatSubList <- c(cellchatSubList, aCellChat)
aCellChat <- readRDS("CellChat_subsetDeligated.rds")
cellchatSubList <- c(cellchatSubList, aCellChat)
aCellchatSubMerge <- readRDS("CellChat_subsetMerge.rds")

levels(cellchatSubList[[1]]@idents)
nlevels(cellchatSubList[[1]]@idents)

sourceCluster <- c(1, 2, 7, 10)
for(i in sourceCluster){
  for(j in 1:nlevels(cellchatSubList[[1]]@idents)){
    netVisual_bubble_MLD(aCellchatSubMerge, source = c(i), target = c(j), apath = "netVisualBubble2_subset")
  }
}

#Generate the files for NATMI on seurat with sub-endotheal cluster IDs, NATMI needs to be run on the commandline
seuratList <- c()
seuratList <- c(seuratList, subset(aSeurat_subset, protocol == "Mock"))
seuratList <- c(seuratList, subset(aSeurat_subset, protocol == "Ligated"))
seuratList <- c(seuratList, subset(aSeurat_subset, protocol == "Deligated"))

for(i in seuratList){
  #From NATMI
  write.csv(100 * (exp(as.matrix(GetAssayData(object = i, assay = "SCT", slot = "data"))) - 1), 
            sprintf("data-%s.csv", i@meta.data$protocol[1]), row.names = T)
  write.csv(Idents(object = i), sprintf("metadata-%s.csv", i@meta.data$protocol[1]), row.names = T)
}

  seuratList <- c()
  seuratList <- c(seuratList, subset(aSeurat_subsetTMP, protocol == "Mock"))
  seuratList <- c(seuratList, subset(aSeurat_subsetTMP, protocol == "Ligated"))
  seuratList <- c(seuratList, subset(aSeurat_subsetTMP, protocol == "Deligated"))
  
  for(i in seuratList){
    #From NATMI
    write.csv(100 * (exp(as.matrix(GetAssayData(object = i, assay = "SCT", slot = "data"))) - 1), 
              sprintf("data2-%s.csv", i@meta.data$protocol[1]), row.names = T)
    write.csv(Idents(object = i), sprintf("meta2data-%s.csv", i@meta.data$protocol[1]), row.names = T)
}

  head(cellchatSubList[[indexList$mock]]@netP$pathways, 10)
  
pathway.union <- union(cellchatSubList[[indexList$mock]]@netP$pathways, 
                       cellchatSubList[[indexList$ligated]]@netP$pathways)
pathway.union <- union(pathway.union, 
                       cellchatSubList[[indexList$deligated]]@netP$pathways)

pathway.union.short <- union(head(cellchatSubList[[indexList$mock]]@netP$pathways, 10), 
                       head(cellchatSubList[[indexList$ligated]]@netP$pathways, 10))
pathway.union.short <- union(pathway.union.short, 
                       head(cellchatSubList[[indexList$deligated]]@netP$pathways, 10))
  
pdf("netAnalysis_signalingRole_heatmap_sub_mock.pdf", height = 9)
  print(netAnalysis_signalingRole_heatmap(cellchatSubList[[indexList$mock]], pattern = "all", 
                                        signaling = pathway.union, title = "Mock", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
pdf("netAnalysis_signalingRole_heatmap_sub_ligated.pdf", height = 9)
  print(netAnalysis_signalingRole_heatmap(cellchatSubList[[indexList$ligated]], pattern = "all", 
                                        signaling = pathway.union, title = "Ligated", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
pdf("netAnalysis_signalingRole_heatmap_sub_deligated.pdf", height = 9)
  print(netAnalysis_signalingRole_heatmap(cellchatSubList[[indexList$deligated]], pattern = "all", 
                                        signaling = pathway.union, title = "Deligated", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()

pdf("netAnalysis_signalingRole_heatmap_sub_short_mock.pdf", height = 9)
print(netAnalysis_signalingRole_heatmap(cellchatSubList[[indexList$mock]], pattern = "all", 
                                        signaling = pathway.union.short, title = "Mock", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
pdf("netAnalysis_signalingRole_heatmap_sub_short_ligated.pdf", height = 9)
print(netAnalysis_signalingRole_heatmap(cellchatSubList[[indexList$ligated]], pattern = "all", 
                                        signaling = pathway.union.short, title = "Ligated", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
pdf("netAnalysis_signalingRole_heatmap_sub_short_deligated.pdf", height = 9)
print(netAnalysis_signalingRole_heatmap(cellchatSubList[[indexList$deligated]], pattern = "all", 
                                        signaling = pathway.union.short, title = "Deligated", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
  


#Reading and subsetting NATMI data---------------------

NATMIFolder <- "JDRSubset2MockvLigated"
fileName <- sprintf("%s/%s/Delta_edges_lrc2p/DOWN-regulated_mean.csv", projectPath, NATMIFolder)
NATMImockvligated <- read.csv(fileName)

fileName <- sprintf("%s/%s/Delta_edges_lrc2p/UP-regulated_mean.csv", projectPath, NATMIFolder)
NATMIligatedvmock <- read.csv(fileName)

NATMIFolder <- "JDRSubset2MockvDeligated"
fileName <- sprintf("%s/%s/Delta_edges_lrc2p/DOWN-regulated_mean.csv", projectPath, NATMIFolder)
NATMImockvdeligated <- read.csv(fileName)

fileName <- sprintf("%s/%s/Delta_edges_lrc2p/UP-regulated_mean.csv", projectPath, NATMIFolder)
NATMIdeligatedvmock <- read.csv(fileName)

NATMIFolder <- "JDRSubset2LigatedvDeligated"
fileName <- sprintf("%s/%s/Delta_edges_lrc2p/DOWN-regulated_mean.csv", projectPath, NATMIFolder)
NATMIligatedvdeligated <- read.csv(fileName)

fileName <- sprintf("%s/%s/Delta_edges_lrc2p/UP-regulated_mean.csv", projectPath, NATMIFolder)
NATMIligatedvdeligated <- read.csv(fileName)

#Up in MockvLigated
ligand <- "Sema7a"
target <- "Fibroblast"
NATMITable <- NATMIligatedvmock
aNATMISubset <- subset(NATMITable, Sending.cluster == "Artery" & Ligand.symbol == ligand & Target.cluster == target)[,c(1, 2, 3, 4, 39)]
aNATMIArtery <- aNATMISubset[order(aNATMISubset$Target.cluster, aNATMISubset$Ligand.symbol),]
aNATMISubset <- subset(NATMITable, Sending.cluster == "Capillary" & Ligand.symbol == ligand & Target.cluster == target)[,c(1, 2, 3, 4, 39)]
aNATMICapillary <- aNATMISubset[order(aNATMISubset$Target.cluster, aNATMISubset$Ligand.symbol),]
aNATMISubset <- subset(NATMITable, Sending.cluster == "Lymphatic" & Ligand.symbol == ligand & Target.cluster == target)[,c(1, 2, 3, 4, 39)]
aNATMILymph <- aNATMISubset[order(aNATMISubset$Target.cluster, aNATMISubset$Ligand.symbol),]
aNATMISubset <- subset(NATMITable, Sending.cluster == "Vein" & Ligand.symbol == ligand & Target.cluster == target)[,c(1, 2, 3, 4, 39)]
aNATMIVein <- aNATMISubset[order(aNATMISubset$Target.cluster, aNATMISubset$Ligand.symbol),]



#Homeostatic-------------------------

homeostaticFile <- "Z:/Next_Generation_Sequencing_Data/scRNAseq/Exp_09202022_Ligated1130_Mock215_Deligated222and308_AA/MergedDatasets/Cnt1_Cnt2_Merged_(SEURAT_v3)_01.rds"
aSeuratHomeostatic <- readRDS(homeostaticFile)

pdf("UMAP_homeostatic.pdf", width = 40, height = 30)
print(DimPlot(aSeuratHomeostatic, reduction = "umap", label = TRUE, label.size = 40, pt.size = 0.5, repel = TRUE) +
        NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()

Idents(aSeuratHomeostatic) <- "seurat_clusters"
levels(aSeuratHomeostatic)
newClusterIDsLong <- c("0_Epithelial", "1_Endothelial", "2_Endothelial", "3_Epithelial", "4_Immune", "5_Immune", 
                       "6_Immune", "7_Endothelial", "8_Epithelial", "9_Endothelial", "10_Fibroblast", 
                       "11_Immune", "12_Immune", "13_Immune", "14_Endothelial", "15_RBC", 
                       "16_Immune", "17_Pericyte", "18_Immune", "19_Immune", "20_Endothelial", 
                       "21_Epithelial", "22_Stroma", "23_Immune", "24_Immune", "25_Endothelial",
                       "26_Fibroblast")
newClusterIDsShort <- c("Epithelial", "Endothelial", "Endothelial", "Epithelial", "Immune", "Immune", 
                        "Immune", "Endothelial", "Epithelial", "Endothelial", "Fibroblast", 
                        "Immune", "Immune", "Immune", "Endothelial", "RBC", 
                        "Immune", "Pericyte", "Immune", "Immune", "Endothelial", 
                        "Epithelial", "Stroma", "Immune", "Immune", "Endothelial",
                        "Fibroblast")

Idents(aSeuratHomeostatic) <- "seurat_clusters"
newClusterIDs <- newClusterIDsShort
names(newClusterIDs) <- levels(aSeuratHomeostatic)
aSeuratHomeostatic <- RenameIdents(aSeuratHomeostatic, newClusterIDs)
levels(aSeuratHomeostatic)

pdf("UMAP_homeostatic_NewID.pdf", width = 40, height = 30)
print(DimPlot(aSeuratHomeostatic, reduction = "umap", label = TRUE, label.size = 40, pt.size = 0.5, repel = TRUE) +
        NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()

kaluckaFile <- "Z:/Amber/2023/02242023_Kalucka2020_ECanalysis/Kalucka_SMG_Integrated_OGRANIDENTITY_(SEURAT_v3)_01.rds"
aSeuratkalucka <- readRDS(kaluckaFile)
levels(aSeuratkalucka)

newClusterIDsShort <- c("Salivary Gland", "Brain", "Heart", "Colon", "Small Intestine", "Kidney", "Liver",
                        "Lung", "EDL Muscle", "Soleus Muscle", "Spleen", "Testis")

newClusterIDs <- newClusterIDsShort
names(newClusterIDs) <- levels(aSeuratkalucka)
aSeuratkalucka <- RenameIdents(aSeuratkalucka, newClusterIDs)
levels(aSeuratkalucka)

pdf("UMAP_kalucka.pdf", width = 40, height = 30)
print(DimPlot(aSeuratkalucka, reduction = "umap", label = FALSE, label.size = 40, pt.size = 0.5, 
              repel = TRUE, order = "Salivary Gland") +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 120)))
dev.off()


