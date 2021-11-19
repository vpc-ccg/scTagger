library(yaml)
library(Seurat)
library(DoubletFinder)

source("functions/preprocessing.R")
source("functions/clustering.R")


args <- commandArgs(trailingOnly = TRUE)
inputPath = args[2]
resultPath = args[3]
config = yaml.load_file(args[1])


preprocessingRes <- preprocessing(inputPath=inputPath, doVlnPlot=config$doVlnPlot,
                                  resultPath=config$resultPath, percentMt=config$percentMt,
                                  nFeatuerRNAMax=config$nFeatuerRNAMax,
                                  nFeatuerRNAMin=config$nFeatuerRNAMin, doSave=config$preDoSave,
                                  project=config$project, nCellMin=config$nCellMin,
                                  vlnPltName=config$vlnPltName, dataName=config$preDataName)

clusteringRes <- clustering(preprocessingRes$data, doEblowPlot=config$doEblowPlot,
                            maxDim=config$maxDim, doUMAP=config$doUMAP, doSave=config$doSave,
                            resultPath=resultPath, nfeatures=config$nfeatures,
                            elbowPlotName=config$elbowPlotName, UMAPName=config$UMAPName,
                            dataName=config$dataName, doDoublet=config$doDoublet)