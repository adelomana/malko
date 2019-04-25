# (This may take a few minutes)
biocLite(c("GEOquery"))
library(GEOquery)
geoFile <- getGEOSuppFiles("GSE60361", makeDirectory=FALSE)
gzFile <- grep("Expression", basename(rownames(geoFile)), value=TRUE)
txtFile <- gsub(".gz", "", gzFile)
gunzip(gzFile, destname=txtFile, remove=TRUE)

library(data.table)
geoData <- fread(txtFile, sep="\t")
geneNames <- unname(unlist(geoData[,1, with=FALSE]))
exprMatrix <- as.matrix(geoData[,-1, with=FALSE])
rm(geoData)
dim(exprMatrix)
rownames(exprMatrix) <- geneNames
exprMatrix[1:5,1:4]

# Remove file downloaded:
file.remove(txtFile)

cellLabels <- paste(file.path(system.file('examples', package='AUCell')), "mouseBrain_cellLabels.tsv", sep="/")
cellLabels <- read.table(cellLabels, row.names=1, header=TRUE, sep="\t")
cellLabels <- as.data.frame(cellLabels)
colnames(cellLabels) <- "CellType"

exprMatrix <- exprMatrix[unique(rownames(exprMatrix)),] # Remove duplicated rows
dim(exprMatrix)

library(SingleCellExperiment)
sceMouseBrain <- SingleCellExperiment(assays = list(counts = exprMatrix),
                                      colData=data.frame(cellLabels[colnames(exprMatrix),, drop=FALSE]))

# setwd("SCENIC_MouseBrain")
dir.create("data")
save(sceMouseBrain, file="data/sceMouseBrain.RData")
