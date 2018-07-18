

# SCENIC Set-up -----------------------------------------------------------

source("https://bioconductor.org/biocLite.R")
biocLite(c("GENIE3", "RcisTarget", "AUCell"))

# Check versions. Correct versions are: AUCell 1.2.4, RcisTarget 1.0.2, and GENIE3 1.2.1 or posterior.
packageVersion("AUCell")
packageVersion("RcisTarget")
packageVersion("GENIE3")

##If you have trouble installing these packages with biocLite, try direct installation from website:
#install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/AUCell_1.2.4.tar.gz", repos=NULL)
#install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/RcisTarget_1.0.2.tar.gz", repos=NULL)
#install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/GENIE3_1.2.1.tar.gz", repos=NULL)

## Install recommended packages for visualization and parallel processing
# Recommended to run AUCell:
biocLite(c("mixtools"))
# To visualize the binary matrices and perform t-SNEs:
biocLite(c("NMF", "Rtsne", "R2HTML"))
# To support paralell execution:
biocLite(c("doMC", "doRNG"))
# To visualize in http://scope.aertslab.org
install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR")

## Install SCENIC

setwd("/Users/MattWall/Desktop/network_package/data/")

# install.packages("devtools")
devtools::install_github("aertslab/SCENIC")


## Download motif databases for RcisTarget (slow: can take >30 min)

# mc9nr: Motif collection version 9: 24k motifs
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")

dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
  descrURL <- gsub(".feather$", ".descr", featherURL)
  if(file.exists(descrURL)) download.file(descrURL, destfile=basename(descrURL))
}

# End SCENIC Set-up -------------------------------------------------------



# SCENIC Analysis ---------------------------------------------------------

# Suppress loading messages when building the HTML
suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(SingleCellExperiment)
})

# To build a personalized report, update this working directory:
install.packages("knitr")
knitr::opts_knit$set(root.dir = 'SCENIC_MouseBrain')

# Create directory to save files generated during run
dir.create("SCENIC_MouseBrain")
setwd("SCENIC_MouseBrain") # Or `knitr::opts_knit$set(root.dir = 'example_results/SCENIC_MouseBrain')` in the first chunk if running a notebook

# Load single-cell dataset


source("https://bioconductor.org/biocLite.R")
biocLite(c("data.table"))
library(data.table)

source("https://bioconductor.org/biocLite.R")
biocLite(c("data.table"))
library(data.table)
setwd("/Users/MattWall/Desktop/network_package/data/")
scData <- fread("malignant.8kgenes.data.csv", sep=",")
geneNames <- unname(unlist(scData[,1, with=FALSE]))
exprMatrix <- as.matrix(scData[,-1, with=FALSE])
rm(scData)
dim(exprMatrix)
rownames(exprMatrix) <- geneNames
exprMatrix[1:5,1:4]

exprMatrix <- exprMatrix[unique(rownames(exprMatrix)),] # Remove duplicated rows
dim(exprMatrix)

# Label cell types
scMetadata <- read.csv("malignant.8kgenes.tumorMetadata.csv", row.names=1,header=TRUE)
colnames(scMetadata) <- "CellType"
cellLabels <- as.data.frame(scMetadata)


# Create single Bioconductor object using SingleCellExperiment
source("https://bioconductor.org/biocLite.R")
biocLite(c("SingleCellExperiment"))

library(SingleCellExperiment)
scHumanMelanoma <- SingleCellExperiment(assays = list(counts = exprMatrix),
                                      colData=data.frame(cellLabels[colnames(exprMatrix),, drop=FALSE]))

# save scHumanMelanoma bioconductor object
dir.create("SCENIC_data")
save(scHumanMelanoma, file="data/scHumanMelanoma.RData")

# Load bioconductor object into Matrix form for SCENIC analysis
load("SCENIC_data/scHumanMelanoma.RData")
exprMat <- counts(scHumanMelanoma)
dim(exprMat)

# Collect information used in plot labels
cellInfo <- colData(scHumanMelanoma)
cellInfo$nGene <- colSums(exprMat>0)
cellInfo <- data.frame(cellInfo)
head(cellInfo)
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=setNames(c("forestgreen", "darkorange", "magenta4", "hotpink", "red3", "skyblue"), 
                                  c("Mel78", "Mel79", "Mel80", "Mel81", "Mel88", "Mel89")))
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

# Set inputs for SCENIC run

library(SCENIC)
org="hgnc" # or hgnc, or dmel
dbDir="/Users/MattWall/Desktop/network_package/data/RcisTarget" # RcisTarget databases location
myDatasetTitle="SCENIC example on Human Melanoma" # choose a name for your analysis
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=4) 

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Databases:
scenicOptions@settings$dbs <- setNames("hg19-tss-centered-10kb-7species.mc9nr.feather", "hg19-500bp-upstream-7species.mc9nr.feather")
scenicOptions@settings$db_mcVersion <- "v9"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# Analyze matrix for noise removal
nCellsPerGene <- apply(exprMat, 1, function(x) sum(x>0))
nCountsPerGene <- apply(exprMat, 1, sum)

summary(nCellsPerGene)
summary(nCountsPerGene)
max(exprMat)
sum(exprMat>0) / sum(exprMat==0)

# Filter 1: Remove cells that on average have fewer than 1 log2(TPM+1) in 1% of cells

minReads <- 2*.01*ncol(exprMat)
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)]
length(genesLeft_minReads)

#Filter 2: Remove genes that appear in fewer than 1% of cells
minSamples <- ncol(exprMat)*.01
nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
length(genesLeft_minCells)

#Filter 3: Only keep the genes in the reference databases

library(RcisTarget)
motifRankings <- importRankings(getDatabases(scenicOptions)[[1]]) # either one, they should have the same genes
genesInDatabase <- colnames(getRanking(motifRankings))
genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)

# Check whether any relevant gene / potential gene of interest is missing:
#interestingGenes <- c("Neurod1", "Sox10", "Dlx1")
#interestingGenes[which(!interestingGenes %in% genesLeft_minCells_inDatabases)]

genesKept <- genesLeft_minCells_inDatabases
saveRDS(genesKept, file=getIntName(scenicOptions, "genesKept"))

exprMat_filtered <- exprMat[genesKept, ]

# To avoid confusion in the coming steps:
rm(exprMat)

## Perform Correlation
corrMat <- cor(t(exprMat_filtered), method="spearman")
# (Only the rows for TFs will be needed needed):
allTFs <- getDbTfs(scenicOptions)
corrMat <- corrMat[which(rownames(corrMat) %in% allTFs),]
saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))



# Run SCENIC --------------------------------------------------------------


# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# Load bioconductor object into Matrix form for SCENIC analysis
# setwd("/Users/MattWall/Desktop/network_package/data/")
load("SCENIC_data/scHumanMelanoma.RData")
exprMat <- counts(scHumanMelanoma)
dim(exprMat)

#Set SCENIC parameters
library(SCENIC)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 4
scenicOptions@settings$seed <- 123

install.packages('zoo')
#Run SCENIC using wrappers
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat)


# End SCENIC run ----------------------------------------------------------


# Visualize results with t-SNE --------------------------------------------

scenicOptions@settings$seed <- 123 # same seed for all of them
# Grid search of different t-SNE settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(20,30,40), perpl=c(20,30,40))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(20,30,40), perpl=c(20,30,40), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot gridsearch t-SNE results as pdf (individual files in int/): 
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE)

par(mfcol=c(3,3))
par(mar=c(1,1,1,1)) #This prevents the plotting error "Error in plot.new() : figure margins too large"
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

# Using only "high-confidence" regulons (normally similar)
par(mfcol=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

# Set tSNE defaults 
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 30
scenicOptions@settings$defaultTsne$perpl <- 40
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


# scenicOptions@settings$devType="png"
runSCENIC_4_aucell_binarize(scenicOptions)

# prepare for plotting
install.packages("rbokeh")
setwd("/Users/MattWall/Desktop/network_package/data/")
logMat <- exprMat # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat) #default t-SNE
savedSelections <- shiny::runApp(aucellApp)

#AUCell_plotTSNE() #to save static plots:
  
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

library(KernSmooth)
library(RColorBrewer)
par(mar=c(1,1,1,1)) #This prevents the plotting error "Error in plot.new() : figure margins too large"
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("XBP1", "ETV1", "TP53")],], plots="Expression")
# Save all AUC into one PDF:
install.packages("Cairo")
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()


# Show several regulons simultaneously

#par(bg = "black")
#par(mfrow=c(1,2))
#regulonNames <- c( "ETV1","XBP1") # replace with the names of the regulons to be viewed
#cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
#text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
#text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)


# Load Regulons -----------------------------------------------------------

regulons <- loadInt(scenicOptions, "regulons")
regulons[c("ETV1", "XBP1")]

# Analyze regulon

install.packages("gplots")
library("gplots")

regulonName <- "MEF2A_extended"
regulonGenes <- regulons[c(regulonName)][[1]]
regulonExpression <- zScore[regulonGenes,]
pdf(paste(paste(regulonName,"heatmap",sep = "_"),"pdf",sep = "."))
heatmap.2(regulonExpression, trace="none")
dev.off()

hist(colVars(etv1Expression))

# Search and display regulon motifs
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
regulonTargetsInfo <- RcisTarget::addLogo(regulonTargetsInfo, motifCol="bestMotif")
regulonTargetsInfo$Genie3Weight <- signif(regulonTargetsInfo$Genie3Weight, 2)

install.packages("DT")
colsToShow <- c("TF", "gene", "nMotifs", "bestMotif", "logo", "NES", "highConfAnnot", "Genie3Weight")
DT::datatable(regulonTargetsInfo[TF=="ETV1" & highConfAnnot==TRUE, colsToShow, with=F], escape=FALSE, filter="top")

######

#check modules from GENIE3

#tfModules_asDF <- loadInt(scenicOptions, "tfModules_asDF")

######

