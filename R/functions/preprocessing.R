preprocessing <- function(inputPath, doVlnPlot=FALSE, resultPath, percentMt=100,
                          nFeatuerRNAMax=2000, nFeatuerRNAMin=200, doSave=TRUE,
                          project="NOV", nCellMin=3,
                          vlnPltName="original_features.jpg",
                          dataName="data.rds"){
  ### Definition: In this function, we read the count matrix and then filter data and normalized it.
  ### Input:
  #### inputPath: Path to the directory of matrix data that cellranger generate it.
  #### doVlnPlot: A Boolean that shows that we need the violent plot or not.
  #### resultPath: Path to the directory that we want to save data there.
  #### percentMt: The threshold that we eliminate the data that have mitochondrial percentage above this number.
  #### nFeatuerRNAMax: Include cells where at last these many features are detected
  #### nFeatuerRNAMin: Include cells where at least this many features are detected
  #### doSave: A boolean number that shows we want to save data in resultPath or not
  #### project: name of the project
  #### nCellMin: min.cells in SeuratObject (Include features detected in at least this many cells. Will subset the counts' matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.)
  #### vlnPltName: Name of violate plot that save in result path
  #### dataName: Name of data that save in result path.
  data.matrix <- Read10X(data.dir = inputPath)

  data <- CreateSeuratObject(counts = data.matrix, project = project,
                             min.cells = nCellMin, min.features = nFeatuerRNAMin)

  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

  if(doVlnPlot){
    features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
    vlnplot <- VlnPlot(data, features = features, ncol = length(features), pt.size=0)
    jpeg(paste0(resultPath, vlnPltName))
    print(vlnplot)
    dev.off()
  }

  data <- subset(data, subset = nFeature_RNA > nFeatuerRNAMin &
                   nFeature_RNA < nFeatuerRNAMax & percent.mt < percentMt)

  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  if(doSave){
    saveRDS(data, file = paste0(resultPath, dataName))
  }
  return (list(data=data, plot=vlnplot))
}
