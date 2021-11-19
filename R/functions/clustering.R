clustering <- function(data, doEblowPlot=FALSE, maxDim=10,
                       doUMAP=TRUE, doSave=TRUE, resultPath, nfeatures=2000,
                       elbowPlotName="elbow.jpg", UMAPName="clustering.jpg",
                       dataName="data-cluster.rds", doDoublet=FALSE,
                       vlnPltName="filter_features.jpg"){
  ### Definition: This function make a clustering plot.
  ### Input:
  #### data: A seurat object.
  #### doEblowPlot: A Boolean that shows that we need the violent plot or not.
  #### maxDim: Dimensions of reduction to use as input
  #### doUMAP:  A Boolean that shows that we need the UMAP plot or not.
  #### doSave: A boolean number that shows we want to save data in resultPath or not
  #### resultPath: Path to the directory that we want to save data there.
  #### nfeatures: Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'
  #### elbowPlotName: Name of elbow plot that save in result path
  #### UMAPName: Name of UMAP plot that save in result path
  #### dataName: Name of data that save in result path.

  data <- FindVariableFeatures(data, selection.method = "vst",
                               nfeatures = nfeatures)
  all.genes <- rownames(data)
  data <- ScaleData(data, features = all.genes)
  data <- RunPCA(data, features = VariableFeatures(object = data))

  if(doEblowPlot){
    elbowPlot <- ElbowPlot(data)

    jpeg(paste0(resultPath, elbowPlotName))
    print(elbowPlot)
    dev.off()
  }

  data <- FindNeighbors(data, dims = 1:maxDim)
  data <- FindClusters(data, resolution = 0.5)
  data <- RunUMAP(data, dims = 1:maxDim)
  print(dim(data))
  if(doDoublet){

    nExp <- round(ncol(data) * (2.3 - 0.8) /200000 * ncol(data))
    print(nExp)
    data <- doubletFinder_v3(data, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:maxDim)
    DF.name = colnames(data@meta.data)[grepl("DF.classification", colnames(data@meta.data))]
    print(nExp)
    jpeg(paste0(resultPath, "Doublet.jpg"))
    doublet <- VlnPlot(data, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
    print(doublet)
    dev.off()
    data = data[, data@meta.data[, DF.name] == "Singlet"]
    Idents(data) <- "orig.ident"
print(colnames(data))
    print(head(Idents(data)))
    features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

    vlnplot <- VlnPlot(data, features = features,
                       ncol = length(features), pt.size=0)
    jpeg(paste0(resultPath, vlnPltName))
    print(vlnplot)
    dev.off()
Idents(data) <- "seurat_clusters"
  }
  print(dim(data))
  if(doUMAP){

    jpeg(paste0(resultPath, UMAPName))
    cluster<-DimPlot(data, reduction = "umap", label=TRUE)
    print(cluster)
    dev.off()
  }

  if(doSave){
    saveRDS(data, file = paste0(resultPath, dataName))
  }
  return(list(data=data))
}