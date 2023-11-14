
if (!require("pacman")) install.packages("pacman")

pacman::p_load(devtools)
devtools::install_github("renozao/NMF@devel")

pacman::p_load(SingleCellExperiment, celda, Seurat, ggplot2, cowplot, scater, patchwork, 
               dplyr, sctransform, glmGamPoi, stringr, scCustomize, scDblFinder, singleCellTK, biomaRt,
               clustree, hdf5r, rhdf5, sctransform, remotes)

install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)

devtools::install_github("sqjin/CellChat")
library("CellChat")

process10x <- function(seuratObject, projName = "project", percentMT = 25, clusterResolution = 1, clustreeP = TRUE){
  
  #debugging assingments
  if(0)  {
    seuratObject = aSeurat_ligated_all
    projName = "AmberLigatedMerge_CH"
    percentMT = 25
    clusterResolution = 1
    clusttreeP = TRUE
  }
  
  dir.create(projName)
  
  seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern = "^mt-")
  
  pdf(paste(projName, "/Feature_Count_MT_Pre.pdf", sep = ""), width = 40, height = 30)
  print(VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) +
          theme(axis.title = element_text(size = 60), axis.text = element_text(size = 60), title = element_text(size = 60))
  )
  dev.off()
  
  seuratObject <- subset(seuratObject, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < percentMT)
  
  pdf(paste(projName, "/Feature_Count_MT_Post.pdf", sep = ""), width = 40, height = 30)
  print(VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) +
          theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40)) 
  ) 
  dev.off()
  
  seuratObject <- NormalizeData(seuratObject)
  seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)
  geneList <- rownames(seuratObject)
  seuratObject <- ScaleData(seuratObject, feature = geneList)
  seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
  
  aplot <-ElbowPlot(seuratObject, ndims = 50)
  pdf(paste(projName, "/ElbowPlot.pdf", sep = ""), width = 40, height = 30)
  print(aplot)
  dev.off()
  
  seuratObject <- FindNeighbors(seuratObject, dims = 1:40)
  
  if(clustreeP){
    seuratTree <- FindClusters(seuratObject, resolution = (seq(0, (clusterResolution*2), by = (clusterResolution/5))))
  
  pdf(paste(projName, "/clusttree.pdf", sep = ""), width = 40, height = 30)
  print(clustree(seuratTree, prefix = "RNA_snn_res."))
  dev.off()
  }
  
  seuratObject <- FindClusters(seuratObject, resolution = clusterResolution)
  
  seuratObject <- RunTSNE(seuratObject, dims = 1:40, check_duplicates=FALSE)
  seuratObject <- RunUMAP(seuratObject, dims = 1:40)
  
  pdf(paste(projName, "/tsne.pdf", sep = ""), width = 40, height = 30)
  print(DimPlot(seuratObject, reduction = "tsne", label = TRUE))
  dev.off()
  
  pdf(paste(projName, "/umap.pdf", sep = ""), width = 40, height = 30)
  print(DimPlot(seuratObject, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
          NoLegend() +
          theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40)))
  dev.off()
  
  featureList = ""
  tryList = ""
  
  baseList <- c("Mucl2", "Aqp5", "Pecam1", "Cdh5", "Prol1", "Sox9", "Sox2", "Sox10", "Krt7", 
               "Krt5", "Krt14", "Trp63", "Acta2", "Kcnip4", "Dcpp1", "Dcpp2", "Dcpp3")
  featureList <- baseList
  
  featureList <- tryCatch(
    {
      sum(GetAssayData(object = seuratObject, slot = "data")["tdTomato",]>0)
      featureList = c("tdTomato", baseList)
    },
    error=function(cond)
    {
      message("Error handler:Dotplot List")
      message(cond)
      tryList = baseList
      return(tryList)
    },
    warning=function(cond)
    {
      message("Warning handler:Dotplot List")
      message(cond)
      tryList = baseList
      return(tryList)
    },
    finally = {}
  )
  
  pdf(paste(projName, "/DotPlot.pdf", sep = ""), width = 40, height = 30)
  
  print(DotPlot(seuratObject, features = featureList, cols = c("lightblue", "red3"), 
                col.min = 0, dot.scale = 20) +
          theme(axis.title = element_text(size = 60), axis.text = element_text(size = 60), 
                legend.text = element_text(size = 60), legend.key.size = unit(3, 'cm'), 
                legend.title  = element_text(size = 60), legend.spacing = unit(3, 'cm')) +
          RotatedAxis())
  dev.off()
  
  return(seuratObject)
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
Read_CellBender_h5_Mat <- function(file_name, use.names = TRUE, unique.features = TRUE) {
  # Check hdf5r installed
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    cli_abort(message = c("Please install hdf5r to read HDF5 files",
                          "i" = "`install.packages('hdf5r')`")
    )
  }
  # Check file
  if (!file.exists(file_name)) {
    stop("File not found")
  }
  
  if (use.names) {
    feature_slot <- 'features/name'
  } else {
    feature_slot <- 'features/id'
  }
  
  # Read file
  infile <- hdf5r::H5File$new(filename = file_name, mode = "r")
  
  counts <- infile[["matrix/data"]]
  indices <- infile[["matrix/indices"]]
  indptr <- infile[["matrix/indptr"]]
  shp <- infile[["matrix/shape"]]
  features <- infile[[paste0("matrix/", feature_slot)]][]
  barcodes <- infile[["matrix/barcodes"]]
  
  
  sparse.mat <- sparseMatrix(
    i = indices[] + 1,
    p = indptr[],
    x = as.numeric(x = counts[]),
    dims = shp[],
    repr = "T"
  )
  
  if (unique.features) {
    features <- make.unique(names = features)
  }
  
  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
  
  infile$close_all()
  
  return(sparse.mat)
}
removeDoublets <- function(seuratObject){
  #debug assignments
  if(0){
    seuratObject <- aSeuratLigated1_filtered
  }
  aSCE <- as.SingleCellExperiment(seuratObject)
  aSCE <- scDblFinder(aSCE, dbr = 0.2)
  
  seuratObject$scDblFinder.score <- aSCE$scDblFinder.score
  seuratObject$scDblFinder.class <- aSCE$scDblFinder.class
  
  dblSeurat <- seuratObject
  dblSeurat <- NormalizeData(dblSeurat)
  dblSeurat <- FindVariableFeatures(dblSeurat, selection.method = "vst", nfeatures = 2000)
  geneList <- rownames(dblSeurat)
  dblSeurat <- ScaleData(dblSeurat, feature = geneList)
  dblSeurat <- RunPCA(dblSeurat, features = VariableFeatures(object = dblSeurat))
  dblSeurat <- RunUMAP(dblSeurat, dims = 1:40)
  
  dblCount <- as.integer(length(rownames(dblSeurat@meta.data))*.2)
  dblSeurat <- doubletFinder_v3(dblSeurat, PCs = 1:40, pN = 0.25, pK = 0.1, nExp = dblCount, reuse.pANN = FALSE, sct = FALSE)
  DF <- sprintf("DF.classifications_0.25_0.1_%i", dblCount)
  DFIndex <- which(colnames(dblSeurat@meta.data) == DF)
  
  seuratObject$DF.classifications <- dblSeurat@meta.data$DF
  
  doublets1 <- subset(seuratObject, scDblFinder.class == "doublet")
  trueDoublets <- subset(doublets1, DF.classifications == "Doublet")
  seuratObject$True_Doublet <- ifelse(colnames(seuratObject)%in% colnames(trueDoublets), "Doublet", "Singlet")
  
  print(sprintf("%i doublets found in %i total cells (%.2f%%)", ncol(trueDoublets), ncol(seuratObject), 
        (ncol(trueDoublets)/ncol(seuratObject)*100)))
  
  seuratObject <- subset(seuratObject, True_Doublet == "Singlet")
  return(seuratObject)
}
loadForIntegrate <- function(fileNamePath, cellBenderP = TRUE, projectName = "Project"){
  #debug assignments
  if(0)
  {
    #fileNamePath <- "Z:/Kenney/E034 Amber reviewer comments/ligated_14_repeat/ligated_14_repeat_CellBender_filtered.h5"
    fileNamePath <- "Z:/Next_Generation_Sequencing_Data/scRNAseq/Exp_022222_5DayDeligated/ML_5_5DayDeligated_022222/outs/filtered_feature_bc_matrix"
    cellBenderP <- FALSE
    projectName <- "Project"
  }
  
  if(cellBenderP){
    aMatrix <- Read_CellBender_h5_Mat(fileNamePath)
    aName <- str_sub(basename(fileNamePath), 1, (str_length(basename(fileNamePath))-12))
    aSeurat <- CreateSeuratObject(aMatrix, project = aName, min.cells = 3, min.genes = 200)
  }else{
    a10X <- Read10X(data.dir = fileNamePath)
    aName <- projectName
    aSeurat <- CreateSeuratObject(a10X, project = aName, min.cells = 3, min.genes = 200)
  }
  aSeurat[["percent.mt"]] <- PercentageFeatureSet(aSeurat, pattern = "^mt-")
  aSeurat <- subset(aSeurat, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  
  aSeurat <- removeDoublets(aSeurat)
  return(aSeurat)
}

#dataPath <- "Z:/Next_Generation_Sequencing_Data/scRNAseq/Exp_021522_Mock_Repeat"
dataPath <- "Z:/Kenney/E034 Amber reviewer comments"
rawDataPath <- "Z:/Next_Generation_Sequencing_Data/scRNAseq"
projectPath <- "Z:/Kenney/E034 Amber reviewer comments"

setwd(projectPath)
set.seed(12354)
options(future.globals.maxSize = 10 * 1024^3)

#Temp dir assign for SCTransformed/Integrated/Cellchat work
workingDir <- "FinalDataSetsSCTransformIntegrated"
dir.create(workingDir)
setwd(sprintf("%s/%s", projectPath, workingDir))

#Cellbender cleanup, all used datasets integrated -------------------------------------
{
  workingDir <- "FinalDataSetsIntegrated"
  dir.create(workingDir)
  setwd(sprintf("%s/%s", projectPath, workingDir))
  
  seuratList <- c()
  aName <- "ligated_14_repeat"
  afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  aSeurat <- loadForIntegrate(afile)
  aSeurat@meta.data[, "protocol"] <- "Ligated"
  seuratList <- c(seuratList, aSeurat)
  
  # aName <- "ligated_14_CH1"
  # afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  # aSeurat <- loadForIntegrate(afile)
  # aSeurat@meta.data[, "protocol"] <- "Ligated"
  # seuratList <- c(seuratList, aSeurat)
  
  # aName <- "ligated_14_CH2"
  # afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  # aSeurat <- loadForIntegrate(afile)
  # aSeurat@meta.data[, "protocol"] <- "Ligated"
  # seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_repeat"
  afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  aSeurat <- loadForIntegrate(afile)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_new1"
  afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  aSeurat <- loadForIntegrate(afile)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_new2"
  afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  aSeurat <- loadForIntegrate(afile)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  # aName <- "mock_14_CH"
  # afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  # aSeurat <- loadForIntegrate(afile)
  # aSeurat@meta.data[, "protocol"] <- "Mock"
  # seuratList <- c(seuratList, aSeurat)
  
  aName <- "deligated"
  afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  aSeurat <- loadForIntegrate(afile)
  aSeurat@meta.data[, "protocol"] <- "Deligated"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "deligated_repeat"
  afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  aSeurat <- loadForIntegrate(afile)
  aSeurat@meta.data[, "protocol"] <- "Deligated"
  seuratList <- c(seuratList, aSeurat)
  
  seuratList <- lapply(X = seuratList, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

  features <- SelectIntegrationFeatures(object.list = seuratList)
  anchors <- FindIntegrationAnchors(object.list = seuratList, anchor.features = features)
  aSeuratIntegrated <- IntegrateData(anchorset = anchors)
  DefaultAssay(aSeuratIntegrated) <- "integrated"
  
  aSeuratIntegrated <- ScaleData(aSeuratIntegrated, verbose = FALSE)
  aSeuratIntegrated <- RunPCA(aSeuratIntegrated, npcs = 40, verbose = FALSE)
  aSeuratIntegrated <- RunUMAP(aSeuratIntegrated, reduction = "pca", dims = 1:40)
  aSeuratIntegrated <- FindNeighbors(aSeuratIntegrated, reduction = "pca", dims = 1:40)
  
  clusterResolution <- 0.6
  seuratTree <- FindClusters(aSeuratIntegrated, resolution = (seq(0, (clusterResolution*2), by = (clusterResolution/5))))
  pdf("All_integrated_clustree.pdf", width = 40, height = 30)
  print(clustree(seuratTree, prefix = "integrated_snn_res.", node_colour = "sc3_stability"))
  dev.off()

  aSeuratIntegrated <- FindClusters(aSeuratIntegrated, resolution = clusterResolution)

  pdf("UMAP_integrate_All.pdf", width = 40, height = 30)
  print(DimPlot(aSeuratIntegrated, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
  dev.off()

  DefaultAssay(aSeuratIntegrated) <- "RNA"
  Idents(aSeuratIntegrated) <- "seurat_clusters"
  aSeuratIntegrated <- NormalizeData(aSeuratIntegrated)
  aSeuratIntegrated <- FindVariableFeatures(aSeuratIntegrated, selection.method = "vst", nfeatures = 2000)
  geneList <- rownames(aSeuratIntegrated)
  aSeuratIntegrated <- ScaleData(aSeuratIntegrated, feature = geneList)
  
  saveRDS(aSeuratIntegrated, file = "AmberAllIntegrated.rds")
  aSeuratIntegrated <- readRDS("AmberAllIntegrated.rds")
  
  saveClusterMarkers(aSeuratIntegrated, projName = "ClusterMarkers")
}

#Cellbender cleanup, all used datasets SCTransform integrated -------------------------------------
{
  workingDir <- "FinalDataSetsSCTransformIntegrated"
  dir.create(workingDir)
  setwd(sprintf("%s/%s", projectPath, workingDir))
  
  seuratList <- c()
  aName <- "ligated_14_repeat"
  afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  aSeurat <- loadForIntegrate(afile)
  aSeurat@meta.data[, "protocol"] <- "Ligated"
  seuratList <- c(seuratList, aSeurat)
  
  # aName <- "ligated_14_CH1"
  # afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  # aSeurat <- loadForIntegrate(afile)
  # aSeurat@meta.data[, "protocol"] <- "Ligated"
  # seuratList <- c(seuratList, aSeurat)
  
  # aName <- "ligated_14_CH2"
  # afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  # aSeurat <- loadForIntegrate(afile)
  # aSeurat@meta.data[, "protocol"] <- "Ligated"
  # seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_repeat"
  afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  aSeurat <- loadForIntegrate(afile)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_new1"
  afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  aSeurat <- loadForIntegrate(afile)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_new2"
  afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  aSeurat <- loadForIntegrate(afile)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  # aName <- "mock_14_CH"
  # afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  # aSeurat <- loadForIntegrate(afile)
  # aSeurat@meta.data[, "protocol"] <- "Mock"
  # seuratList <- c(seuratList, aSeurat)
  
  aName <- "deligated"
  afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  aSeurat <- loadForIntegrate(afile)
  aSeurat@meta.data[, "protocol"] <- "Deligated"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "deligated_repeat"
  afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
  aSeurat <- loadForIntegrate(afile)
  aSeurat@meta.data[, "protocol"] <- "Deligated"
  seuratList <- c(seuratList, aSeurat)
  
  seuratList <- lapply(X = seuratList, FUN = function(x) {
    x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
  })
  
  features <- SelectIntegrationFeatures(object.list = seuratList, nfeatures = 3000)
  seuratList <- PrepSCTIntegration(object.list = seuratList, anchor.features = features)
  anchors <- FindIntegrationAnchors(object.list = seuratList, anchor.features = features, normalization.method = "SCT")
  aSeuratIntegrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  DefaultAssay(aSeuratIntegrated) <- "integrated"
  
  aSeuratIntegrated <- ScaleData(aSeuratIntegrated, verbose = FALSE)
  aSeuratIntegrated <- RunPCA(aSeuratIntegrated, npcs = 40, verbose = FALSE)
  aSeuratIntegrated <- RunUMAP(aSeuratIntegrated, reduction = "pca", dims = 1:40)
  aSeuratIntegrated <- FindNeighbors(aSeuratIntegrated, reduction = "pca", dims = 1:40)
  
  clusterResolution <- 0.6
  seuratTree <- FindClusters(aSeuratIntegrated, resolution = (seq(0, (clusterResolution*2), by = (clusterResolution/5))))
  pdf("All_integrated_clustree.pdf", width = 40, height = 30)
  print(clustree(seuratTree, prefix = "integrated_snn_res.", node_colour = "sc3_stability"))
  dev.off()
  
  aSeuratIntegrated <- FindClusters(aSeuratIntegrated, resolution = clusterResolution)
  
  Idents(aSeuratIntegrated) <- "seurat_clusters"
  levels(aSeuratIntegrated)
  pdf("UMAP_integrate_All.pdf", width = 40, height = 30)
  print(DimPlot(aSeuratIntegrated, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
          #NoLegend() +
          theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
  dev.off()
  
  Idents(aSeuratIntegrated) <- "protocol"
  levels(aSeuratIntegrated)
  pdf("UMAP_integrate_All_protocol.pdf", width = 40, height = 30)
  print(DimPlot(aSeuratIntegrated, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
          #NoLegend() +
          theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
  dev.off()
  
  # DefaultAssay(aSeuratIntegrated) <- "RNA"
  # Idents(aSeuratIntegrated) <- "seurat_clusters"
  # aSeuratIntegrated <- NormalizeData(aSeuratIntegrated)
  # aSeuratIntegrated <- FindVariableFeatures(aSeuratIntegrated, selection.method = "vst", nfeatures = 2000)
  # geneList <- rownames(aSeuratIntegrated)
  # aSeuratIntegrated <- ScaleData(aSeuratIntegrated, feature = geneList)
  
  DefaultAssay(aSeuratIntegrated) <- "SCT"
  Idents(aSeuratIntegrated) <- "seurat_clusters"
  aSeuratIntegrated <- PrepSCTFindMarkers(aSeuratIntegrated, assay = "SCT", verbose = TRUE)
  
  saveRDS(aSeuratIntegrated, file = "AmberAllSCTransformIntegrated.rds")
  aSeuratIntegrated <- readRDS("AmberAllSCTransformIntegrated.rds")
  
  saveClusterMarkers(aSeuratIntegrated, projName = "ClusterMarkers")
  
  levels(aSeuratIntegrated)
  newClusterIDsLong <- c("0_NK", "1_Endothelial", "2_Endothelial", "3_Endothelial", "4_Fibroblast", "5_TCell", 
                         "6_Endothelial", "7_Duct", "8_End_Immune", "9_NK", "10_Macrophage", "11_Endothelial", 
                         "12_Endothelial", "13_Smooth_Muscle", "14_Smooth_Muscle", "15_Endothelial", "16_Glial", "17_TCell",
                         "18_TCell", "19", "20_Immune", "21_Immune", "22_Acinar", "23_BCell", "24", 
                         "25")
  newClusterIDsShort <- c("NK", "Endothelial", "Endothelial", "Endothelial", "Fibroblast", "TCell", 
                          "Endothelial", "Duct", "End_Immune", "NK", "Macrophage", "Endothelial", 
                          "Endothelial", "Smooth_Muscle", "Smooth_Muscle", "Endothelial", "Glial", "TCell",
                          "TCell", "19", "Immune", "Immune", "Acinar", "BCell", "24", 
                          "25")
  Idents(aSeuratIntegrated) <- "seurat_clusters"
  newClusterIDs <- newClusterIDsShort
  names(newClusterIDs) <- levels(aSeuratIntegrated)
  aSeuratIntegrated <- RenameIdents(aSeuratIntegrated, newClusterIDs)
  levels(aSeuratIntegrated)
  
  aseuratLigated <- subset(aSeuratIntegrated, protocol == "Ligated")
  aseuratMock <- subset(aSeuratIntegrated, protocol == "Mock")
  aseuratDeligated <- subset(aSeuratIntegrated, protocol == "Deligated")
  
  #Generate the files for NATMI, NATMI needs to be run on the commandline
  seuratList <- c()
  seuratList <- c(seuratList, subset(aSeuratIntegrated, protocol == "Ligated"))
  seuratList <- c(seuratList, subset(aSeuratIntegrated, protocol == "Mock"))
  seuratList <- c(seuratList, subset(aSeuratIntegrated, protocol == "Deligated"))

  for(i in seuratList){
  
  #From Web
  # write.table(as.matrix(GetAssayData(object = i, slot = "counts")), 
  #             sprintf("counts-%s.csv", i@meta.data$protocol[1]), 
  #             sep = ',', row.names = T, col.names = T, quote = F)
    # write.table(as.matrix(i@active.ident), sprintf('annotation-%s.csv', i@meta.data$protocol[1]), 
    #             sep = ',', row.names = T, col.names = T, quote = F)
    
  #From NATMI
    write.csv(100 * (exp(as.matrix(GetAssayData(object = i, assay = "SCT", slot = "data"))) - 1), 
              sprintf("data-%s.csv", i@meta.data$protocol[1]), row.names = T)
    write.csv(Idents(object = i), sprintf("metadata-%s.csv", i@meta.data$protocol[1]), row.names = T)
    
  }
  
  #Cellchat
  data.input <- GetAssayData(aseuratDeligated, assay = "SCT", slot = "data")
  labels <- Idents(aseuratDeligated)
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
  
  #groupSize <- as.numeric(table(acellchat@idents))
  
  acellchat@meta$protocol <- "Deligated"
  saveRDS(acellchat, file = "AmberDeligatedSCTransformCellChat.rds")
  acellchat <- readRDS("AmberDeligatedSCTransformCellChat.rds")
  
  acellchat@meta$protocol <- "Ligated"
  saveRDS(acellchat, file = "AmberLigatedSCTransformCellChat.rds")
  acellchat <- readRDS("AmberLigatedSCTransformCellChat.rds")
  
  acellchat@meta$protocol <- "Mock"
  saveRDS(acellchat, file = "AmberMockSCTransformCellChat.rds")
  acellchat <- readRDS("AmberMockSCTransformCellChat.rds")
  
  cellchatList <- c()
  acellchat <- readRDS("AmberDeligatedSCTransformCellChat.rds")
  cellchatList <- c(cellchatList, acellchat)
  
  acellchat <- readRDS("AmberLigatedSCTransformCellChat.rds")
  cellchatList <- c(cellchatList, acellchat)
        
  acellchat <- readRDS("AmberMockSCTransformCellChat.rds")
  cellchatList <- c(cellchatList, acellchat)
  
}
#Cellbender Integrated dotplot and Cellchat ---------------------------
{
  workingDir <- "All_Integrated"
  dir.create(workingDir)
  
  aSeuratIntegrated <- readRDS("C:/Dropbox/Postdoc/TempStorage/AmberAllIntegrated.rds")
  aSeuratIntegrated <- readRDS("Z:/Kenney/E034 Amber reviewer comments/All_Integrated/AmberAllIntegrated.rds")
  
  levels(aSeuratIntegrated)
  newClusterIDsLong <- c("0_Endothelial", "1_NK", "2_Endothelial", "3_Epithelial", "4_Fibroblast", "5_Endothelial", 
                         "6_Macrophage", "7_Endothelial", "8_TCell", "9_Epithelial", "10_NK", "11_Endothelial", 
                         "12_Epithelial", "13_TCell", "14_Endothelial", "15_Macrophage", "16_Pericyte", "17_Smooth Muscle",
                         "18_Endothelial", "19_RBC", "20_Epithelial", "21_Macrophage", "22_Immune", "23_Pericyte", "24", 
                         "25_Macrophage", "26")
  newClusterIDsShort <- c("Endothelial", "NK", "Endothelial", "Epithelial", "Fibroblast", "Endothelial", 
                          "Macrophage", "Endothelial", "TCell", "Epithelial", "NK", "Endothelial", 
                          "Epithelial", "TCell", "Endothelial", "Macrophage", "Pericyte", "Smooth Muscle",
                          "Endothelial", "RBC", "Epithelial", "Macrophage", "Immune", "Pericyte", "24", 
                          "Macrophage", "26")
  Idents(aSeuratIntegrated) <- "seurat_clusters"
  newClusterIDs <- newClusterIDsShort
  names(newClusterIDs) <- levels(aSeuratIntegrated)
  aSeuratIntegrated <- RenameIdents(aSeuratIntegrated, newClusterIDs)
  levels(aSeuratIntegrated)
  
  pdf(sprintf("./%s/UMAP_integrate_All_newID.pdf", workingDir), width = 40, height = 30)
  print(DimPlot(aSeuratIntegrated, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
          #NoLegend() +
          theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
  dev.off()
  
  VlnPlot(aSeuratIntegrated, features = c("Pdgfra", "Pdgfrb"))
  VlnPlot(aSeuratIntegrated, features = c("Pdgfrb", "Acta2"))
  VlnPlot(aSeuratIntegrated, features = c("Col1a1"))
  VlnPlot(aSeuratIntegrated, features = c("Ptprc"))
  VlnPlot(aSeuratIntegrated, features = c("Acta2"))
  VlnPlot(aSeuratIntegrated, features = c("Pecam1"))
  
  DefaultAssay(aSeuratIntegrated) <- "RNA"
  aSeuratFibroblast <- subset(aSeuratIntegrated, idents = "Fibroblast")
  aSeuratEndothelial <- subset(aSeuratIntegrated, idents = "Endothelial")
  
  Idents(aSeuratIntegrated) <- "orig.ident"
  Idents(aSeuratFibroblast) <- "orig.ident"
  Idents(aSeuratEndothelial) <- "orig.ident"
  
  for(i in levels(aSeuratIntegrated)){
    aSeuratFib <- subset(aSeuratFibroblast, idents = i)
    aSeurat <- subset(aSeuratIntegrated, idents = i)
    print(sprintf("%i cells in %s of %i total", ncol(aSeuratFib), i, ncol(aSeurat)))
  }
  
  Idents(aSeuratFibroblast) <- "protocol"
  levels(aSeuratFibroblast)
  DefaultAssay(aSeuratIntegrated) <- "RNA"
  DefaultAssay(aSeuratEndothelial) <- "RNA"
  DefaultAssay(aSeuratFibroblast) <- "RNA"
  
  featureList <- c("Abi3bp", "Aebp1", "Bgn", "Cilp", "Col11a1", "Col12a1", "Col16a1", "Col28a1", "Col4a3", "Col4a4", 
                   "Col4a5", "Col8a1", "Col8a2", "Creld2", "Efemp2", "Eln", "Emilin1", "Fbn2", "Fmod", "Gas6",
                   "Igfbp7", "Igfbp6", "Ltbp2", "Mfap2", "Mfap4", "Mfap5", "Mgp", "Omd", "Postn", "Slit2", "Sparc", "Spp1",
                   "Svep1", "Thbs2", "Thbs4", "Acta2", "Col1a1", "Dcn", "Dpp4", "Fap", "Fn1", "Pdgfra", "Pdgfrb", "Vim",
                   "Gli1", "Gli3", "Smo")
  
  endFeatureList <- c("Pecam1", "Cdh5", "Kdr", "Flt1", "Adgre4", "Eng", "Emcn", "Col4a1", "Timp3", "Plvap", "Igfbp7", 
                     "Ly6c1", "Fabp4", "Egfl7", "Cxcl12", "Sparc", "Ly6a", "Id3", "Cav1", "Esam")
  
  
  aSeuratFibroblastTmp <- aSeuratFibroblast
  aSeuratFibroblastTmp@active.ident <- factor(aSeuratFibroblastTmp@active.ident, 
      levels = c("0_Ligated", "0_Mock", "0_Deligated", "1_Ligated", "1_Mock", "1_Deligated",
                 "3_Ligated", "3_Mock", "3_Deligated", "4_Ligated", "4_Mock", "4_Deligated", "2_Mock"))
  levels(aSeuratFibroblastTmp)
  
  pdf(sprintf("./%s/Dotplot Endothelial Fig 6 reorder.pdf", workingDir), width = 40, height = 30)
  print(DotPlot(aSeuratEndothelial, features = endFeatureList, cols = c("lightblue", "red3"), 
                col.min = 0, dot.scale = 20) +
          theme(axis.title = element_text(size = 60), axis.text = element_text(size = 30), 
                legend.text = element_text(size = 60), legend.key.size = unit(3, 'cm'), 
                legend.title  = element_text(size = 60), legend.spacing = unit(3, 'cm')) +
          RotatedAxis())
  dev.off()
  
  
  aSeuratFibroblast <- NormalizeData(aSeuratFibroblast)
  aSeuratFibroblast <- FindVariableFeatures(aSeuratFibroblast, selection.method = "vst", nfeatures = 2000)
  aSeuratFibroblast <- ScaleData(aSeuratFibroblast, verbose = FALSE)
  aSeuratFibroblast <- RunPCA(aSeuratFibroblast, npcs = 40, verbose = FALSE)
  aSeuratFibroblast <- RunUMAP(aSeuratFibroblast, reduction = "pca", dims = 1:40)
  aSeuratFibroblast <- FindNeighbors(aSeuratFibroblast, reduction = "pca", dims = 1:40)
  
  clusterResolution <- 0.5
  seuratTree <- FindClusters(aSeuratFibroblast, resolution = (seq(0, (clusterResolution*2), by = (clusterResolution/5))))
  pdf(sprintf("./%s/clustree Fibroblast.pdf", workingDir), width = 40, height = 30)
  print(clustree(seuratTree, prefix = "RNA_snn_res."))
  dev.off()
  
  aSeuratFibroblast <- FindClusters(aSeuratFibroblast, resolution = 0.1)
  aSeuratFibroblast$subCluster <- paste0(aSeuratFibroblast$seurat_clusters, "_", aSeuratFibroblast$protocol)
  Idents(aSeuratFibroblast) <- "subCluster"
  
  VlnPlot(aSeuratFibroblast, features = c("Pdgfra", "Pdgfrb"))
  
  aSeuratEndothelial <- NormalizeData(aSeuratEndothelial)
  aSeuratEndothelial <- FindVariableFeatures(aSeuratEndothelial, selection.method = "vst", nfeatures = 2000)
  aSeuratEndothelial <- ScaleData(aSeuratEndothelial, verbose = FALSE)
  aSeuratEndothelial <- RunPCA(aSeuratEndothelial, npcs = 40, verbose = FALSE)
  aSeuratEndothelial <- RunUMAP(aSeuratEndothelial, reduction = "pca", dims = 1:40)
  aSeuratEndothelial <- FindNeighbors(aSeuratEndothelial, reduction = "pca", dims = 1:40)
  
  pdf(sprintf("./%s/endothelial heatmap.pdf", workingDir), width = 40, height = 30)
  print(DoHeatmap(aSeuratEndothelial,  features = VariableFeatures(aSeuratEndothelial), group.by = "orig.ident"))
  dev.off()
 
  Idents(aSeuratIntegrated) <- "orig.ident"
  Idents(aSeuratEndothelial) <- "orig.ident"
  levels(aSeuratEndothelial)
  
  aSeuratEndothelial <- ScaleData(object = aSeuratEndothelial, features = rownames(aSeuratEndothelial))
  
  for(i in levels(aSeuratEndothelial)){
  print(i)
  allMarkers <- FindMarkers(aSeuratEndothelial, ident.1 = i, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, 
                               group.by = 'orig.ident')
  
  allMarkers <- arrange(allMarkers, desc(abs(allMarkers[,2])))
  topMarkers <- allMarkers[1:100,]
  
  pdf(sprintf("./%s/endothelial %s heatmap.pdf", workingDir, i), width = 40, height = 30)
  print(DoHeatmap(aSeuratEndothelial,  features = rownames(topMarkers), group.by = "orig.ident")) 
  }
  

  
}

#No Cellbender cleanup, all datasets integrated --------------------
{
  workingDir <- "All_Integrated_Raw"
  dir.create(workingDir)
  
  seuratList <- c()
  aName <- "ligated_14_repeat"
  localDataPath <- "Exp_216_Ligated_Repeat_02/ML_1_ligated_113021"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Ligated"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "ligated_14_CH1"
  localDataPath <- "Exp_231_LigatedMock_noFGF2organoid/Ligated/ML-1"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Ligated"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "ligated_14_CH2"
  localDataPath <- "Exp_231_LigatedMock_noFGF2organoid/Ligated/ML-4"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Ligated"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_repeat"
  localDataPath <- "Exp_021522_Mock_Repeat/ML_2_Mock_021522"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_new1"
  localDataPath <- "Exp_06132023_MockSurgeryReplicate/Mock0613_1"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_new2"
  localDataPath <- "Exp_06132023_MockSurgeryReplicate/Mock0613_2"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_CH"
  localDataPath <- "Exp_231_LigatedMock_noFGF2organoid/Mock/ML-3"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "deligated"
  localDataPath <- "Exp_022222_5DayDeligated/ML_5_5DayDeligated_022222"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Deligated"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "deligated_repeat"
  localDataPath <- "Exp_030822_5DayDeligatedRepeat/ML-6_5DayDeligated_030822"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Deligated"
  seuratList <- c(seuratList, aSeurat)
  
  seuratList <- lapply(X = seuratList, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  features <- SelectIntegrationFeatures(object.list = seuratList)
  anchors <- FindIntegrationAnchors(object.list = seuratList, anchor.features = features)
  aSeuratIntegrated <- IntegrateData(anchorset = anchors)
  DefaultAssay(aSeuratIntegrated) <- "integrated"
  
  aSeuratIntegrated <- ScaleData(aSeuratIntegrated, verbose = FALSE)
  aSeuratIntegrated <- RunPCA(aSeuratIntegrated, npcs = 40, verbose = FALSE)
  aSeuratIntegrated <- RunUMAP(aSeuratIntegrated, reduction = "pca", dims = 1:40)
  aSeuratIntegrated <- FindNeighbors(aSeuratIntegrated, reduction = "pca", dims = 1:40)
  
  clusterResolution <- 1
  seuratTree <- FindClusters(aSeuratIntegrated, resolution = (seq(0, (clusterResolution*2), by = (clusterResolution/5))))
  pdf(sprintf("./%s/All_integrated_clustree.pdf", workingDir), width = 40, height = 30)
  print(clustree(seuratTree, prefix = "integrated_snn_res.", node_colour = "sc3_stability"))
  dev.off()
  
  aSeuratIntegrated <- FindClusters(aSeuratIntegrated, resolution = 0.6)
  
  saveRDS(aSeuratIntegrated, file = "./%s/AmberAllIntegrated.rds", workingDir)
  #E033_decontX_all <- readRDS("./All_Integrated/AmberAllIntegrated.rds")
  
  tmp <- subset(aSeuratIntegrated, orig.ident == "mock_14_new2_CellBender")
  
  pdf(sprintf("./%s/UMAP_integrate_All.pdf", workingDir), width = 40, height = 30)
  print(DimPlot(aSeuratIntegrated, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
          #NoLegend() +
          theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
  dev.off()
  
  Idents(aSeuratIntegrated) <- "orig.ident"
  for(i in levels(aSeuratIntegrated)){
    print(i)
    subSeurat <- subset(aSeuratIntegrated, orig.ident == i)
    
    pdf(sprintf("./%s/UMAP_integrate_%s.pdf", workingDir, i), width = 40, height = 30)
    print(DimPlot(subSeurat, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
            #NoLegend() +
            theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
    dev.off()
  }
  Idents(aSeuratIntegrated) <- "seurat_clusters"
  
  saveClusterMarkers(aSeuratIntegrated, projName = "All_Integrated")
}

#No Cellbender Integrated dotplot and Cellchat ----------------
{
  workingDir <- "All_Integrated_Raw"
  dir.create(workingDir)
  
  #saveRDS(aSeuratIntegratedRaw, file = "./%s/AmberAllIntegrated.rds", workingDir)
  aSeuratIntegratedRaw <- readRDS("C:/Dropbox/Postdoc/TempStorage/AmberAllIntegratedraw.rds")
  
  levels(aSeuratIntegratedRaw)
  newClusterIDsLong <- c("0_Endothelial)", "1_Immune", "2_Endothelial", "3_Epithelial", "4_Endothelial", "5_Fibroblast", 
                         "6_Immune", "7_Immune", "8_Endothelial", "9_Epithelial", "10_Immune", "11_Immune", "12_Epithelial",
                         "13_Immune", "14_Endothelial", "15_Smooth Muscle", "16_Smooth Muscle", "17_Endothelial", 
                         "18_Epithelial", "19", "20_Immune", "21_Glial", "22", "23", "24_Immune", "25_Immune", 
                         "26_Immune", "27_Immune", "28_Fibroblast")
  newClusterIDsShort <- c("Endothelial", "Immune", "Endothelial", "Epithelial", "Endothelial", "Fibroblast", 
                          "Immune", "Immune", "Endothelial", "Epithelial", "Immune", "Immune", "Epithelial",
                          "Immune", "Endothelial", "Smooth Muscle", "Smooth Muscle", "Endothelial", 
                          "Epithelial", "19", "Immune", "Glial", "22", "23", "Immune", "Immune", 
                          "Immune", "Immune", "Fibroblast")
  
  Idents(aSeuratIntegratedRaw) <- "seurat_clusters"
  newClusterIDs <- newClusterIDsShort
  names(newClusterIDs) <- levels(aSeuratIntegratedRaw)
  aSeuratIntegratedRaw <- RenameIdents(aSeuratIntegratedRaw, newClusterIDs)
  levels(aSeuratIntegratedRaw)
  
  pdf(sprintf("./%s/UMAP_integrate_All_newIDShort.pdf", workingDir), width = 40, height = 30)
  print(DimPlot(aSeuratIntegratedRaw, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
          #NoLegend() +
          theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
  dev.off()
  
  aSeuratFibroblastRaw <- subset(aSeuratIntegratedRaw, idents = "Fibroblast")
  Idents(aSeuratFibroblastRaw) <- "orig.ident"
  Idents(aSeuratIntegratedRaw) <- "orig.ident"
  levels(aSeuratFibroblastRaw)
  DefaultAssay(aSeuratFibroblastRaw) <- "RNA"
  
  for(i in levels(aSeuratIntegratedRaw)){
    aSeuratFib <- subset(aSeuratFibroblastRaw, idents = i)
    aSeurat <- subset(aSeuratIntegratedRaw, idents = i)
    print(sprintf("%i cells in %s of %i total", ncol(aSeuratFib), i, ncol(aSeurat)))
  }
  
  featureList <- c("Abi3bp", "Aebp1", "Bgn", "Cilp", "Col11a1", "Col12a1", "Col16a1", "Col28a1", "Col4a3", "Col4a4", 
                   "Col4a5", "Col8a1", "Col8a2", "Creld2", "Efemp2", "Eln", "Emilin1", "Fbn2", "Fmod", "Gas6",
                   "Igfbp7", "Igfbp6", "Ltbp2", "Mfap2", "Mfap4", "Mfap5", "Mgp", "Omd", "Postn", "Slit2", "Sparc", "Spp1",
                   "Svep1", "Thbs2", "Thbs4", "Acta2", "Col1a1", "Dcn", "Dpp4", "Fap", "Fn1", "Pdgfra", "Pdgfrb", "Vim",
                   "Gli1", "Gli3", "Smo")
  
  aSeuratFibroblastTmp <- aSeuratFibroblastRaw
  aSeuratFibroblastTmp@active.ident <- factor(aSeuratFibroblastTmp@active.ident, 
                                              levels = c("0_Ligated", "0_Mock", "0_Deligated", "1_Ligated", "1_Mock", "1_Deligated",
                                                         "3_Ligated", "3_Mock", "3_Deligated", "4_Ligated", "4_Mock", "4_Deligated", "2_Mock"))
  levels(aSeuratFibroblastTmp)
  
  pdf(sprintf("./%s/Dotplot Fibroblast Fig 6 origIdent.pdf", workingDir), width = 40, height = 30)
  print(DotPlot(aSeuratFibroblastTmp, features = featureList, cols = c("lightblue", "red3"), 
                col.min = 0, dot.scale = 20) +
          theme(axis.title = element_text(size = 60), axis.text = element_text(size = 30), 
                legend.text = element_text(size = 60), legend.key.size = unit(3, 'cm'), 
                legend.title  = element_text(size = 60), legend.spacing = unit(3, 'cm')) +
          RotatedAxis())
  dev.off()
  
  
  aSeuratFibroblast <- NormalizeData(aSeuratFibroblast)
  aSeuratFibroblast <- FindVariableFeatures(aSeuratFibroblast, selection.method = "vst", nfeatures = 2000)
  aSeuratFibroblast <- ScaleData(aSeuratFibroblast, verbose = FALSE)
  aSeuratFibroblast <- RunPCA(aSeuratFibroblast, npcs = 40, verbose = FALSE)
  aSeuratFibroblast <- RunUMAP(aSeuratFibroblast, reduction = "pca", dims = 1:40)
  aSeuratFibroblast <- FindNeighbors(aSeuratFibroblast, reduction = "pca", dims = 1:40)
  
  clusterResolution <- 0.5
  seuratTree <- FindClusters(aSeuratFibroblast, resolution = (seq(0, (clusterResolution*2), by = (clusterResolution/5))))
  pdf(sprintf("./%s/clustree Fibroblast.pdf", workingDir), width = 40, height = 30)
  print(clustree(seuratTree, prefix = "RNA_snn_res."))
  dev.off()
  
  aSeuratFibroblast <- FindClusters(aSeuratFibroblast, resolution = 0.1)
  aSeuratFibroblast$subCluster <- paste0(aSeuratFibroblast$seurat_clusters, "_", aSeuratFibroblast$protocol)
  Idents(aSeuratFibroblast) <- "subCluster"
  
  VlnPlot(aSeuratFibroblast, features = c("Pdgfra", "Pdgfrb"))
}

#No cellbender Mock and control integrated -------------------
{
  workingDir <- "Ctrl_Mock_Integrated_Raw"
  dir.create(workingDir)
  
  seuratList <- c()
  
  aName <- "mock_14_repeat"
  localDataPath <- "Exp_021522_Mock_Repeat/ML_2_Mock_021522"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_new1"
  localDataPath <- "Exp_06132023_MockSurgeryReplicate/Mock0613_1"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_new2"
  localDataPath <- "Exp_06132023_MockSurgeryReplicate/Mock0613_2"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "mock_14_CH"
  localDataPath <- "Exp_231_LigatedMock_noFGF2organoid/Mock/ML-3"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Mock"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "CtrlStroma_1"
  localDataPath <- "Exp_08082022_C57ControlStromaPrep/ML_12_CtrlStroma_081522"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Control"
  seuratList <- c(seuratList, aSeurat)
  
  aName <- "CtrlStroma_2"
  localDataPath <- "Exp_08082022_C57ControlStromaPrep/ML_13_CtrlStroma_081522"
  afile <- sprintf("%s/%s/outs/filtered_feature_bc_matrix", rawDataPath, localDataPath)
  aSeurat <- loadForIntegrate(afile, cellBenderP = FALSE, projectName = aName)
  aSeurat@meta.data[, "protocol"] <- "Control"
  seuratList <- c(seuratList, aSeurat)
  
  seuratList <- lapply(X = seuratList, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  features <- SelectIntegrationFeatures(object.list = seuratList)
  anchors <- FindIntegrationAnchors(object.list = seuratList, anchor.features = features)
  aSeuratIntegrated <- IntegrateData(anchorset = anchors)
  DefaultAssay(aSeuratIntegrated) <- "integrated"
  
  aSeuratIntegrated <- ScaleData(aSeuratIntegrated, verbose = FALSE)
  aSeuratIntegrated <- RunPCA(aSeuratIntegrated, npcs = 40, verbose = FALSE)
  aSeuratIntegrated <- RunUMAP(aSeuratIntegrated, reduction = "pca", dims = 1:40)
  aSeuratIntegrated <- FindNeighbors(aSeuratIntegrated, reduction = "pca", dims = 1:40)
  
  clusterResolution <- 0.6
  seuratTree <- FindClusters(aSeuratIntegrated, resolution = (seq(0, (clusterResolution*2), by = (clusterResolution/5))))
  pdf(sprintf("./%s/All_integrated_clustree.pdf", workingDir), width = 40, height = 30)
  print(clustree(seuratTree, prefix = "integrated_snn_res.", node_colour = "sc3_stability"))
  dev.off()
  
  aSeuratIntegrated <- FindClusters(aSeuratIntegrated, resolution = clusterResolution)
  
  saveRDS(aSeuratIntegrated, file = sprintf("./%s/Ctrl_Mock_Integrated.rds", workingDir))
  aSeuratIntegrated <- readRDS(sprintf("./%s/Ctrl_Mock_Integrated.rds", workingDir))
  
  pdf(sprintf("./%s/UMAP_integrate_All.pdf", workingDir), width = 40, height = 30)
  print(DimPlot(aSeuratIntegrated, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
          #NoLegend() +
          theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
  dev.off()
  
  Idents(aSeuratIntegrated) <- "orig.ident"
  for(i in levels(aSeuratIntegrated)){
    print(i)
    subSeurat <- subset(aSeuratIntegrated, orig.ident == i)
    
    pdf(sprintf("./%s/UMAP_integrate_%s.pdf", workingDir, i), width = 40, height = 30)
    print(DimPlot(subSeurat, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
            #NoLegend() +
            theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
    dev.off()
  }
  
  DefaultAssay(aSeuratIntegrated) <- "RNA"
  Idents(aSeuratIntegrated) <- "seurat_clusters"
  aSeuratIntegrated <- NormalizeData(aSeuratIntegrated)
  aSeuratIntegrated <- FindVariableFeatures(aSeuratIntegrated, selection.method = "vst", nfeatures = 2000)
  geneList <- rownames(aSeuratIntegrated)
  aSeuratIntegrated <- ScaleData(aSeuratIntegrated, feature = geneList)
  # aSeuratIntegrated <- RunPCA(aSeuratIntegrated, features = VariableFeatures(object = aSeuratIntegrated))
  # aSeuratIntegrated <- FindNeighbors(aSeuratIntegrated, dims = 1:40)
  
  saveRDS(aSeuratIntegrated, file = sprintf("./%s/Ctrl_Mock_Integrated2.rds", workingDir))
  aSeuratIntegrated <- readRDS(sprintf("./%s/Ctrl_Mock_Integrated2.rds", workingDir))
  
  #heatmap of top DE genes
  saveClusterMarkers(aSeuratIntegrated, projName = workingDir)
  
  Idents(aSeuratIntegrated) <- "orig.ident"
  aSeuratEndothelial <- aSeuratIntegrated
  endFeatureList <- c("Pecam1", "Cdh5", "Kdr", "Flt1", "Adgre4", "Eng", "Emcn", "Col4a1", "Timp3", "Plvap", "Igfbp7", 
                      "Ly6c1", "Fabp4", "Egfl7", "Cxcl12", "Sparc", "Ly6a", "Id3", "Cav1", "Esam")
  
  levels(aSeuratEndothelial)
  allMarkers <- FindMarkers(aSeuratEndothelial, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,
                               ident.1 = 'mock_14_repeat', group.by = 'orig.ident')
  allMarkers <- arrange(allMarkers, desc(abs(allMarkers[,2])))
  topMarkers <- allMarkers[1:100,]
  topMarkers <- arrange(topMarkers, desc((topMarkers[,2])))
  print(DoHeatmap(aSeuratEndothelial,  features = rownames(topMarkers), group.by = "orig.ident"))
  
  topMarkers2 <- topMarkers
  markerList <- rbind(topMarkers, topMarkers2)
  markerList2 <- rbind(markerList2, markerList)
  
  markerList = matrix(data = NA, nrow = 0, ncol = 5)
  colnames(markerList) <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
  
  for(i in levels(aSeuratEndothelial)){
    print(i)
    allMarkers <- FindMarkers(aSeuratEndothelial, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                   ident.1 = i, group.by = "orig.ident")
#    allMarkers <- arrange(allMarkers, desc(abs(allMarkers[,2])))
    allMarkers <- arrange(allMarkers, desc(abs(allMarkers[,2])))
    topMarkers <- allMarkers[1:100,]
    topMarkers <- arrange(topMarkers, desc((topMarkers[,2])))
    markerList <- rbind(markerList, topMarkers)
    
    pdf(sprintf("./%s/Heatmap topMarkers %s.pdf", workingDir, i), width = 40, height = 30)
    print(DoHeatmap(aSeuratEndothelial,  features = rownames(topMarkers), group.by = "orig.ident"))
    dev.off()
  }
  
  pdf(sprintf("./%s/Heatmap AllMarkers.pdf", workingDir), width = 40, height = 30)
  print(DoHeatmap(aSeuratEndothelial,  features = rownames(markerList), group.by = "orig.ident"))
  dev.off()
  
  #Get gene list of control EC to compare to mock EC
  levels(aSeuratIntegrated)
  newClusterIDsLong <- c("0_Endothelial", "1_Epithelial", "2_Epithelial", "3_Endothelial", "4_NK", "5_Macrophage", 
                         "6_StriatedDuct", "7_Endothelial", "8_Endothelial", "9_Endothelial", "10_Fibroblast",
                         "11_SMG_Acinar", "12_NK", "13_Pericyte", "14_Duct", "15_Macrophage", "16_Endothelial",
                         "17_TCell", "18", "19_Immune", "20_NK", "21_Epithelial", "22_Macrophage", "23_SMG_Duct",
                         "24_Dividing", "25_Glial", "26_Dividing", "27_TCell", "28", "29_Myoepithelial", "30",
                         "31_Endothelial", "32_BCell", "33")
  newClusterIDsShort <- c("Endothelial", "Epithelial", "Epithelial", "Endothelial", "NK", "Macrophage", 
                          "StriatedDuct", "Endothelial", "Endothelial", "Endothelial", "Fibroblast",
                          "SMG_Acinar", "NK", "Pericyte", "Duct", "Macrophage", "Endothelial",
                          "TCell", "18", "Immune", "NK", "Epithelial", "Macrophage", "SMG_Duct",
                          "Dividing", "Glial", "Dividing", "TCell", "28", "Myoepithelial", "30",
                          "Endothelial", "BCell", "33")
  Idents(aSeuratIntegrated) <- "seurat_clusters"
  newClusterIDs <- newClusterIDsShort
  names(newClusterIDs) <- levels(aSeuratIntegrated)
  aSeuratIntegrated <- RenameIdents(aSeuratIntegrated, newClusterIDs)
  levels(aSeuratIntegrated)

  
  # doublets1 <- subset(seuratObject, scDblFinder.class == "doublet")
  # trueDoublets <- subset(doublets1, DF.classifications == "Doublet")
  # seuratObject$True_Doublet <- ifelse(colnames(seuratObject)%in% colnames(trueDoublets), "Doublet", "Singlet")
  
  aSeuratIntegrated$type_ID <- paste0(aSeuratIntegrated$protocol, "_", aSeuratIntegrated@active.ident)
  aSeuratNoMockEC <- subset(aSeuratIntegrated, type_ID != "Mock_Endothelial")
  
  allMarkers <- FindMarkers(aSeuratNoMockEC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                            ident.1 = "Control_Endothelial", group.by = "type_ID")
  rm(aSeuratNoMockEC) #do not need anymore, free up memory
  
  allMarkers %>% top_n(n = 500, wt = avg_log2FC) -> top500
  
  aSeuratEC <- subset(aSeuratIntegrated, idents = "Endothelial")
  pdf(sprintf("./%s/Heatmap Endothelial Control EC Markers.pdf", workingDir), width = 40, height = 30)
  print(DoHeatmap(aSeuratIntegrated,  features = rownames(top500), group.by = "orig.ident"))
  dev.off()
  
  allMarkers %>% top_n(n = 20, wt = avg_log2FC) -> top20
  
  Idents(aSeuratEC) <- "orig.ident"
  levels(aSeuratEC)
  pdf(sprintf("./%s/Dotplot Endothelial 20 Control EC Markers.pdf", workingDir), width = 40, height = 30)
  print(DotPlot(aSeuratEC, features = rownames(top20), cols = c("lightblue", "red3"), 
                col.min = 0, dot.scale = 20) +
          theme(axis.title = element_text(size = 60), axis.text = element_text(size = 60), 
                legend.text = element_text(size = 60), legend.key.size = unit(3, 'cm'), 
                legend.title  = element_text(size = 60), legend.spacing = unit(3, 'cm')) +
          RotatedAxis())
  dev.off()
  
  Stress_Genes_of_Interest <- c("Atf3", "Jun", "Btg2", "Cdkn1a", "Bad", "Bax")
  pdf(sprintf("./%s/Dotplot Endothelial stress Markers.pdf", workingDir), width = 40, height = 30)
  print(DotPlot(aSeuratEC, features = Stress_Genes_of_Interest, cols = c("lightblue", "red3"), 
                col.min = 0, dot.scale = 20) +
          theme(axis.title = element_text(size = 60), axis.text = element_text(size = 60), 
                legend.text = element_text(size = 60), legend.key.size = unit(3, 'cm'), 
                legend.title  = element_text(size = 60), legend.spacing = unit(3, 'cm')) +
          RotatedAxis())
  dev.off()
}

#Merge Trial ---------------
{
#Not to be used, cells to not cluster 'properly'
# aName <- "ligated_14_old"
# afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
# aMatrix <- Read_CellBender_h5_Mat(afile)
# aSeurat1 <- CreateSeuratObject(aMatrix, project = sprintf("%s_CellBender", aName), min.cells = 3, min.genes = 200)

aName <- "ligated_14_repeat"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aMatrix <- Read_CellBender_h5_Mat(afile)
aSeuratLigated1 <- CreateSeuratObject(aMatrix, project = sprintf("%s_CellBender", aName), min.cells = 3, min.genes = 200)

aName <- "ligated_14_CH1"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aMatrix <- Read_CellBender_h5_Mat(afile)
aSeuratLigated2 <- CreateSeuratObject(aMatrix, project = sprintf("%s_CellBender", aName), min.cells = 3, min.genes = 200)

aName <- "ligated_14_CH2"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aMatrix <- Read_CellBender_h5_Mat(afile)
aSeuratLigated3 <- CreateSeuratObject(aMatrix, project = sprintf("%s_CellBender", aName), min.cells = 3, min.genes = 200)

aSeurat_ligated_all <- merge(aSeuratLigated1, y = c(aSeuratLigated2, aSeuratLigated3), add.cell.ids = c("Repeat", "CH1", "CH2"), 
                            project = "Ligated")

aSeurat_ligated_all <- process10x(aSeurat_ligated_all, projName = "AmberLigatedMerge_CH", percentMT = 25, clusterResolution = 1)

saveRDS(aSeurat_all, file = "AmberLigatedMerge.rds")
aSeurat_all <- readRDS("AmberLigatedMerge.rds")

Idents(aSeurat_ligated_all) <- "orig.ident"
pdf("UMAP_AmberLigatedMerge_CH_origIdent.pdf", width = 40, height = 30)
print(DimPlot(aSeurat_ligated_all, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40)))
dev.off()

aName <- "mock_14_new1"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aMatrix <- Read_CellBender_h5_Mat(afile)
aSeuratMock1 <- CreateSeuratObject(aMatrix, project = sprintf("%s_CellBender", aName), min.cells = 3, min.genes = 200)
    
aName <- "mock_14_new2"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aMatrix <- Read_CellBender_h5_Mat(afile)
aSeuratMock2 <- CreateSeuratObject(aMatrix, project = sprintf("%s_CellBender", aName), min.cells = 3, min.genes = 200)

aName <- "mock_14_repeat"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aMatrix <- Read_CellBender_h5_Mat(afile)
aSeuratMock3 <- CreateSeuratObject(aMatrix, project = sprintf("%s_CellBender", aName), min.cells = 3, min.genes = 200)

aName <- "mock_14_CH"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aMatrix <- Read_CellBender_h5_Mat(afile)
aSeuratMock4 <- CreateSeuratObject(aMatrix, project = sprintf("%s_CellBender", aName), min.cells = 3, min.genes = 200)

#Not to be used, cells to not cluster 'properly'
# aName <- "mock_14_old"
# afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
# aMatrix <- Read_CellBender_h5_Mat(afile)
# aSeurat4 <- CreateSeuratObject(aMatrix, project = sprintf("%s_CellBender", aName), min.cells = 3, min.genes = 200)
}

#Merge workflow ---------------
{
aSeurat_all <- merge(aSeurat1, y = c(aSeurat2, aSeurat3, aSeurat4), add.cell.ids = c("new1", "new2", "old", "repeat"), 
                     project = "Mock")

aSeurat_all <- process10x(aSeurat_all, projName = "AmberMockMerge", percentMT = 25, clusterResolution = 1)

saveRDS(aSeurat_all, file = "AmberMockMerge.rds")
#E033_decontX_all <- readRDS("AmberMockMerge.rds")

Idents(aSeurat_all) <- "orig.ident"
levels(aSeurat_all)
levels(aSeurat_all) <- c("mock_14_repeat_CellBender", "mock_14_old_CellBender", "mock_14_new2_CellBender", "mock_14_new1_CellBender")
pdf("temp_reverse.pdf", width = 40, height = 30)
print(DimPlot(aSeurat_all, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40)))
dev.off()
}


#Integrate Workflow --------------------
{
#intList <- c(aSeurat1, aSeurat2, aSeurat3, aSeurat4)
#intList <- c(aSeurat1, aSeurat2, aSeurat3)
#intList <- c(aSeuratLigated1, aSeuratMock1, aSeuratMock2, aSeuratMock3)

  aSeuratLigated1[["percent.mt"]] <- PercentageFeatureSet(aSeuratLigated1, pattern = "^mt-")
  aSeuratLigated1_filtered <- subset(aSeuratLigated1, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
 
  aSeuratLigated2[["percent.mt"]] <- PercentageFeatureSet(aSeuratLigated2, pattern = "^mt-")
  aSeuratLigated2_filtered <- subset(aSeuratLigated2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  
  aSeuratLigated3[["percent.mt"]] <- PercentageFeatureSet(aSeuratLigated3, pattern = "^mt-")
  aSeuratLigated3_filtered <- subset(aSeuratLigated3, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  
intList <- c(aSeuratLigated1_filtered, aSeuratLigated2_filtered, aSeuratLigated3_filtered)

aSeuratMock1[["percent.mt"]] <- PercentageFeatureSet(aSeuratMock1, pattern = "^mt-")
aSeuratMock1_filtered <- subset(aSeuratMock1, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)

aSeuratMock2[["percent.mt"]] <- PercentageFeatureSet(aSeuratMock2, pattern = "^mt-")
aSeuratMock2_filtered <- subset(aSeuratMock2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)

aSeuratMock3[["percent.mt"]] <- PercentageFeatureSet(aSeuratMock3, pattern = "^mt-")
aSeuratMock3_filtered <- subset(aSeuratMock3, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)

aSeuratMock4[["percent.mt"]] <- PercentageFeatureSet(aSeuratMock4, pattern = "^mt-")
aSeuratMock4_filtered <- subset(aSeuratMock4, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)

intList <- c(aSeuratMock1_filtered, aSeuratMock2_filtered, aSeuratMock3_filtered, aSeuratMock4_filtered)

intList <- lapply(X = intList, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = intList)
anchors <- FindIntegrationAnchors(object.list = intList, anchor.features = features)
aSeuratIntegrated <- IntegrateData(anchorset = anchors)
DefaultAssay(aSeuratIntegrated) <- "integrated"

# saveRDS(aSeuratIntegrated, file = "AmberMockIntegrate.rds")
# aSeuratIntegrated <- readRDS("AmberMockIntegrate.rds")

aSeuratIntegrated <- ScaleData(aSeuratIntegrated, verbose = FALSE)
aSeuratIntegrated <- RunPCA(aSeuratIntegrated, npcs = 40, verbose = FALSE)
aSeuratIntegrated <- RunUMAP(aSeuratIntegrated, reduction = "pca", dims = 1:40)
aSeuratIntegrated <- FindNeighbors(aSeuratIntegrated, reduction = "pca", dims = 1:40)

clusterResolution <- 1
seuratTree <- FindClusters(aSeuratIntegrated, resolution = (seq(0, (clusterResolution*2), by = (clusterResolution/5))))
pdf("mock_integrate_CH_clusttree_stablity.pdf", width = 40, height = 30)
  print(clustree(seuratTree, prefix = "integrated_snn_res.", node_colour = "sc3_stability"))
dev.off()

aSeuratIntegrated <- FindClusters(aSeuratIntegrated, resolution = 0.2)

Idents(aSeuratIntegrated) <- "orig.ident"
Idents(aSeuratIntegrated) <- "seurat_clusters"
levels(aSeuratIntegrated)
levels(aSeuratIntegrated) <- c("mock_14_new1_CellBender", "mock_14_new2_CellBender", "mock_14_old_CellBender", "mock_14_repeat_CellBender")
levels(aSeuratIntegrated) <- c("mock_14_repeat_CellBender", "mock_14_old_CellBender", "mock_14_new2_CellBender", "mock_14_new1_CellBender")

tmp <- subset(aSeuratIntegrated, orig.ident == "mock_14_new2_CellBender")

pdf("UMAP_integrate_mock_CH_new2.pdf", width = 40, height = 30)
print(DimPlot(tmp, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()


pdf("UMAP_mock_14_repeat_CellBender.pdf", width = 40, height = 30)
print(DimPlot(tmp, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()
}

#Mock and ligated mixed -------------------
{
  aSeuratIntegrated@assays$RNA <- aSeuratIntegrated@assays$integrated
  DefaultAssay(aSeuratIntegrated) <- "RNA"
  
aSeurat_all <- merge(aSeuratIntegrated, y = c(aSeuratLigated), add.cell.ids = c("Mock", "Ligated"), 
                     project = "AmberData")

aSeurat_all <- process10x(aSeurat_all, projName = "AmberIntegrateMerge", percentMT = 25, clusterResolution = 1, clustreeP = FALSE)
}

#Cell chat demo -------------------------

apath <- sprintf("%s/ML_2_Mock_021522/outs/filtered_feature_bc_matrix", dataPath)
mock_1 <- Read10X(data.dir = apath)
mock_1 <- CreateSeuratObject(mock_1, project = "mock_1", min.cells = 3, min.features = 200)
mock_1 <- process10x(mock_1, projName = "mock_1")

new_cluster_ids = paste0("C", levels(mock_1))
names(new_cluster_ids) <- levels(mock_1)
mock_1 <- RenameIdents(mock_1, new_cluster_ids)

data.input <- GetAssayData(mock_1, assay = "RNA", slot = "data")
labels <- Idents(mock_1)
meta <- data.frame(group = labels, row.names = names(labels))

acellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

acellchat <- addMeta(acellchat, meta = meta) 
acellchat <- setIdent(acellchat, ident.use = "group") # set "group" as default cell identity
levels(acellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(acellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

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

groupSize <- as.numeric(table(acellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(acellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(acellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")

mat <- acellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("TGFb") 
# Hierarchy plot
avertex.receiver = seq(1,10) # a numeric vector. Does not seem to behave different with different values
par(mfrow=c(1,1))
netVisual_aggregate(acellchat, signaling = pathways.show,  vertex.receiver = avertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(acellchat, signaling = pathways.show, layout = "circle")

netVisual_aggregate(acellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(acellchat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(acellchat, signaling = pathways.show)

pairLR.TGFb <- extractEnrichedLR(acellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(2) # a numeric vector
vertex.receiver <- c()
vertex.receiver <- c(vertex.receiver, 2)
vertex.receiver = c(5, 11, 14, 15) # a numeric vector
netVisual_individual(acellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver, 
                     layout = "hierarchy")
netVisual_individual(acellchat, signaling = pathways.show,  pairLR.use = LR.show, sources.use = c(2))
netVisual_individual(acellchat, signaling = pathways.show,  pairLR.use = LR.show, targets.use = c(2))
#> [[1]]
# Circle plot
netVisual_individual(acellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Chord diagram
netVisual_individual(acellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Access all the signaling pathways showing significant communications
pathways.show.all <- acellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(acellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(acellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(acellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(acellchat, sources.use = c(2), targets.use = c(3), remove.isolate = FALSE)
netVisual_bubble(acellchat, sources.use = c(2), remove.isolate = FALSE)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
print(pathways.show.all)
netVisual_bubble(acellchat, sources.use = 2, targets.use = c(1:15), signaling = c("CXCL","CDH5", "LAMININ"), remove.isolate = FALSE)

# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(acellchat, signaling = c("CXCL","CDH5", "LAMININ"))
netVisual_bubble(acellchat, sources.use = c(2), targets.use = c(1:15), pairLR.use = pairLR.use, remove.isolate = TRUE)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(acellchat, sources.use = 2, targets.use = c(1:15), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions received by C7
netVisual_chord_gene(acellchat, sources.use = c(1:15), targets.use = 2, legend.pos.x = 15)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(acellchat)

netAnalysis_signalingRole_heatmap(acellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(acellchat, pattern = "incoming")
netAnalysis_signalingRole_heatmap(acellchat, pattern = "all")

netAnalysis_river(acellchat, pattern = "outgoing")
netAnalysis_river(acellchat, pattern = "incoming")

netAnalysis_dot(acellchat, pattern = "outgoing")
netAnalysis_dot(acellchat, pattern = "incoming")

for(i in cellchatList){
  pdf(sprintf("./cellchat/netAnalysis_river_outgoing_%s.pdf", i@meta$protocol[1]), width = 40, height = 30)
  print(netAnalysis_river(i, pattern = "outgoing", font.size = 8, font.size.title = 30))
  dev.off()
  pdf(sprintf("./cellchat/netAnalysis_river_incoming_%s.pdf", i@meta$protocol[1]), width = 40, height = 30)
  print(netAnalysis_river(i, pattern = "incoming", font.size = 8, font.size.title = 30))
  dev.off()
  
  pdf(sprintf("./cellchat/netAnalysis_dot_outgoing_%s.pdf", i@meta$protocol[1]), width = 40, height = 30)
  print(netAnalysis_dot(i, pattern = "outgoing", font.size = 20, font.size.title = 30, dot.size = c(1,10)))
  dev.off()
  pdf(sprintf("./cellchat/netAnalysis_dot_incoming_%s.pdf", i@meta$protocol[1]), width = 40, height = 30)
  print(netAnalysis_dot(i, pattern = "incoming", font.size = 20, font.size.title = 30, dot.size = c(1,10)))
  dev.off() 
}

#save dotplots
for(i in cellchatList){
  pdf(sprintf("./cellchat/netAnalysis_dot_outgoing_Endothelial_%s.pdf", i@meta$protocol[1]), width = 40, height = 30)
  print(netAnalysis_dot(i, pattern = "outgoing", font.size = 20, font.size.title = 30, 
                        dot.size = c(1,10), group.show = c("Endothelial")))
  dev.off()
  pdf(sprintf("./cellchat/netAnalysis_dot_incoming_Endothelial_%s.pdf", i@meta$protocol[1]), width = 40, height = 30)
  print(netAnalysis_dot(i, pattern = "incoming", font.size = 20, font.size.title = 30, 
                        dot.size = c(1,10), group.show = c("Endothelial")))
  dev.off() 
}

#save heatmaps
for(i in cellchatList){
  pdf(sprintf("./cellchat/heatmap_count_%s.pdf", i@meta$protocol[1]), width = 40, height = 30)
  print(netVisual_heatmap(i, color.heatmap = "Reds", slot.name = "netP", measure = "count", 
        font.size = 40, font.size.title = 50))
  dev.off()
  
  pdf(sprintf("./cellchat/heatmap_weight_%s.pdf", i@meta$protocol[1]), width = 40, height = 30)
  print(netVisual_heatmap(i, color.heatmap = "Reds", slot.name = "netP", measure = "weight",
        font.size = 40, font.size.title = 50))
  dev.off()
}

#Incoming vs outgoing strength
for(i in cellchatList){
  pdf(sprintf("./cellchat/signalingRole_scatter_%s.pdf", i@meta$protocol[1]))
  print(netAnalysis_signalingRole_scatter(i))
  dev.off()
}

aCellchatMerge <- mergeCellChat(cellchatList, add.names = c(cellchatList[[1]]@meta$protocol[1], 
                                          cellchatList[[2]]@meta$protocol[1], 
                                          cellchatList[[3]]@meta$protocol[1]))

#1 Deligated
#2 Ligated
#3 Mock
deligatedIndex <- 1
ligatedIndex <- 2
mockIndex <- 3
compareInteractions(aCellchatMerge, show.legend = F, group = c(1, 2, 3))
compareInteractions(aCellchatMerge, show.legend = F, group = c(1, 2, 3), measure = "weight")

netVisual_diffInteraction(aCellchatMerge, weight.scale = T, measure = "weight", comparison = c(1, 2))
netVisual_diffInteraction(aCellchatMerge, weight.scale = T, measure = "weight", comparison = c(1, 3))
netVisual_diffInteraction(aCellchatMerge, weight.scale = T, measure = "weight", comparison = c(2, 3))

pdf(sprintf("./cellchat/netVisual_heatmap_MockvLigated.pdf"))
print(netVisual_heatmap(aCellchatMerge, comparison = c(mockIndex, ligatedIndex),  measure = "weight",
                  title.name = "Mock vs Ligated differential interaction strength"))
dev.off()

pdf(sprintf("./cellchat/netVisual_heatmap_MockvDeligated.pdf"))
netVisual_heatmap(aCellchatMerge, comparison = c(mockIndex, deligatedIndex), measure = "weight",
                  title.name = "Mock vs Deligated differential interaction strength")
dev.off()

pdf(sprintf("./cellchat/netVisual_heatmap_LigatedvDeligated.pdf"))
netVisual_heatmap(aCellchatMerge, comparison = c(ligatedIndex, deligatedIndex), measure = "weight",
                  title.name = "Ligated vs Deligated differential interaction strength")
dev.off()

pdf(sprintf("./cellchat/signalingChanges_scatter_MockvLigated.pdf"))
print(netAnalysis_signalingChanges_scatter(aCellchatMerge, idents.use = "Endothelial", 
                                     comparison = c(mockIndex, ligatedIndex)))
dev.off()
pdf(sprintf("./cellchat/signalingChanges_scatter_MockvDeligated.pdf"))
print(netAnalysis_signalingChanges_scatter(aCellchatMerge, idents.use = "Endothelial", 
                                           comparison = c(mockIndex, deligatedIndex)))
dev.off()
pdf(sprintf("./cellchat/signalingChanges_scatter_LigatedvDeligated.pdf"))
print(netAnalysis_signalingChanges_scatter(aCellchatMerge, idents.use = "Endothelial", 
                                           comparison = c(ligatedIndex, deligatedIndex)))
dev.off()

aCellchatMerge <- computeNetSimilarityPairwise(aCellchatMerge, type = "functional")

saveRDS(aCellchatMerge, file = "AmberDeligatedSCTransformCellChatMerge.rds")
aCellchatMerge <- readRDS("AmberDeligatedSCTransformCellChatMerge.rds")

# library(reticulate)
# path_to_python <- "C:/anaconda3"
# use_python(path_to_python)

# reticulate::install_python(version = "3.7:latest")
# reticulate::py_config()
# virtualenv_create("CellChat-env", version = "3.7.9")
# use_virtualenv("CellChat-env")
# reticulate::py_install(packages = "umap-learn")
# reticulate::py_install(packages = 'umap-learn')
# virtualenv_install("CellChat-env", packages = "umap-learn")
# py_list_packages(envname = "CellChat-env")
# tmp <- py_list_packages()
# umap-learn <- import("umap-learn")
# scipy <- import("scipy")

#Identify signaling groups based on functional similarity
aCellchatMerge <- netEmbedding(aCellchatMerge, type = "functional", umap.method = "uwot")
aCellchatMerge <- netClustering(aCellchatMerge, type = "functional", do.parallel = FALSE)
netVisual_embeddingPairwise(aCellchatMerge, type = "functional", label.size = 3.5)

#Identify signaling groups based on structure similarity
aCellchatMerge <- computeNetSimilarityPairwise(aCellchatMerge, type = "structural")
aCellchatMerge <- netEmbedding(aCellchatMerge, type = "structural", umap.method = "uwot")
aCellchatMerge <- netClustering(aCellchatMerge, type = "structural", do.parallel = FALSE)
netVisual_embeddingPairwise(aCellchatMerge, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(aCellchatMerge, type = "structural", nCol = 2)


pdf(sprintf("./cellchat/rankSimilarity_MockvLigated.pdf"))
  print(rankSimilarity(aCellchatMerge, type = "functional", comparison2 = c(mockIndex, ligatedIndex)))
dev.off()

pdf(sprintf("./cellchat/rankSimilarity_MockvDeligated.pdf"))
  rankSimilarity(aCellchatMerge, type = "functional", comparison2 = c(mockIndex, deligatedIndex))
dev.off()

pdf(sprintf("./cellchat/rankSimilarity_LigatedvDeligated.pdf"))
  rankSimilarity(aCellchatMerge, type = "functional", comparison2 = c(ligatedIndex, deligatedIndex))
dev.off()

adata <- netVisual_bubble(aCellchatMerge, sources.use = 2, targets.use = c(3, 4, 5, 7),  
                 comparison = c(mockIndex, ligatedIndex), angle.x = 45, return.data = TRUE)

netVisual_bubble(aCellchatMerge, sources.use = 2, targets.use = c(3, 4, 5, 7),  comparison = c(mockIndex, ligatedIndex), 
                 max.dataset = ligatedIndex, title.name = "Increased signaling in Ligated vs Mock", angle.x = 45, 
                 remove.isolate = T)

#End to Fib communication, Up in first
pdf(sprintf("./cellchat/netVisual_bubble_MockvLigated_End-FB.pdf"))
print(netVisual_bubble(aCellchatMerge, sources.use = 2, targets.use = c(3),  comparison = c(mockIndex, ligatedIndex), 
                 max.dataset = mockIndex, title.name = "Endothelial to Fibroblast Increased signaling in Mock vs Ligated", angle.x = 45, 
                 remove.isolate = T))
dev.off()

pdf(sprintf("./cellchat/netVisual_bubble_LigatedvDeligated_End-FB.pdf"))
print(netVisual_bubble(aCellchatMerge, sources.use = 2, targets.use = c(3),  comparison = c(ligatedIndex, deligatedIndex), 
                       max.dataset = ligatedIndex, title.name = "Endothelial to Fibroblast Increased signaling in Ligated vs Deligated", angle.x = 45, 
                       remove.isolate = T))
dev.off()

pdf(sprintf("./cellchat/netVisual_bubble_MockvDeligated_End-FB.pdf"))
print(netVisual_bubble(aCellchatMerge, sources.use = 2, targets.use = c(3),  comparison = c(mockIndex, deligatedIndex), 
                       max.dataset = mockIndex, title.name = "Endothelial to Fibroblast Increased signaling in Mock vs deligated", angle.x = 45, 
                       remove.isolate = T))
dev.off()

#End to Fib communication, Up in second
pdf(sprintf("./cellchat/netVisual_bubble_LigatedvMock_End-FB.pdf"))
print(netVisual_bubble(aCellchatMerge, sources.use = 2, targets.use = c(3),  comparison = c(mockIndex, ligatedIndex), 
                       max.dataset = ligatedIndex, title.name = "Endothelial to Fibroblast Increased signaling in Ligated vs Mock", angle.x = 45, 
                       remove.isolate = T))
dev.off()

pdf(sprintf("./cellchat/netVisual_bubble_DeligatedvLigated_End-FB.pdf"))
print(netVisual_bubble(aCellchatMerge, sources.use = 2, targets.use = c(3),  comparison = c(ligatedIndex, deligatedIndex), 
                       max.dataset = deligatedIndex, title.name = "Endothelial to Fibroblast Increased signaling in Deligated vs Ligated", angle.x = 45, 
                       remove.isolate = T))
dev.off()

pdf(sprintf("./cellchat/netVisual_bubble_DeligatedvMock_End-FB.pdf"))
print(netVisual_bubble(aCellchatMerge, sources.use = 2, targets.use = c(3, 1, 4),  comparison = c(mockIndex, deligatedIndex), 
                       max.dataset = deligatedIndex, title.name = "Endothelial to Fibroblast Increased signaling in Deligated vs Mock", angle.x = 45, 
                       remove.isolate = T))
dev.off()


tmp <- sort(pathway.union)

# data defines ----------------
# 14 day Mock 215 = 
#   "/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_021522_Mock_Repeat/ML_2_Mock_021522/outs/filtered_feature_bc_matrix"
# 
# 14 day Ligated1130="/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/Exp_216_Ligated_Repeat_02/ML_1_ligated_113021/outs/filtered_feature_bc_matrix"
# 
# New 14 day Mock 1 = "Z:/Next_Generation_Sequencing_Data/scRNAseq/Exp_06132023_MockSurgeryReplicate/Mock0613_1/outs/filtered_feature_bc_matrix"
# 
# New 14 day Mock 2 ="Z:/Next_Generation_Sequencing_Data/scRNAseq/Exp_06132023_MockSurgeryReplicate/Mock0613_2/outs/filtered_feature_bc_matrix"

#   Below are the datasets that were processed with singulomics. Cells were frozen and shipped to Singulomics for Library prep and sequencing.
# 
# Old 14 Day Mock ="/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/March2021_DuctalLigationvsMock/68.195.245.254/results/Mock/filtered_feature_bc_matrix"
# 
# Old 14 Day Ligated="/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/scRNAseq/March2021_DuctalLigationvsMock/68.195.245.254/results/Ligated/filtered_feature_bc_matrix"


#Possible other data
#Z:\Next_Generation_Sequencing_Data\scRNAseq\Exp_231_LigatedMock_noFGF2organoid


