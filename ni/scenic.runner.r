# 1. SCENIC Set-up ----

# 0.1. SCENIC dependencies
source("https://bioconductor.org/biocLite.R")
#biocLite(c("GENIE3", "RcisTarget", "AUCell"))

library("GENIE3")
library("RcisTarget")
library("AUCell")

# Check versions. Correct versions are: AUCell 1.2.4, RcisTarget 1.0.2, and GENIE3 1.2.1 or posterior.
packageVersion("AUCell")
packageVersion("RcisTarget")
packageVersion("GENIE3")

## Install recommended packages for visualization and parallel processing
# Recommended to run AUCell:
#biocLite(c("mixtools"))
library("mixtools")
# To visualize the binary matrices and perform t-SNEs:
#biocLite(c("NMF", "Rtsne", "R2HTML"))
library("NMF")
library("Rtsne")
library("R2HTML")
# To support paralell execution:
#biocLite(c("doMC", "doRNG"))
library("doMC")
library("doRNG")
# To visualize in http://scope.aertslab.org
#install.packages("devtools")
library(devtools)
devtools::install_github("aertslab/SCopeLoomR")

library("RCurl")
library("XML")

# 0.2. Install SCENIC
setwd("/Volumes/omics4tb/alomana/projects/mscni/src/scenic/")

# install.packages("devtools")
devtools::install_github("aertslab/SCENIC")

# 0.3. Download databases
## Download motif databases for RcisTarget (slow: can take >30 min)

dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs

dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
  descrURL <- gsub(".feather$", ".descr", featherURL)
  if(file.exists(descrURL)) download.file(descrURL, destfile=basename(descrURL))
}

# End SCENIC Set-up

# 1. download test data

# (This may take a few minutes)
#biocLite(c("GEOquery"))
library(GEOquery)
geoFile <- getGEOSuppFiles("GSE60361", makeDirectory=FALSE)
gzFile <- grep("Expression", basename(rownames(geoFile)), value=TRUE)
txtFile <- gsub(".gz", "", gzFile)
gunzip(gzFile, destname=txtFile, remove=TRUE)

library(data.table)
geoData <- fread(txtFile, sep="\t")
musGeneNames <- unname(unlist(geoData[,1, with=FALSE]))
musExprMatrix <- as.matrix(geoData[,-1, with=FALSE])
rm(geoData)
dim(musExprMatrix)
rownames(musExprMatrix) <- musGeneNames
musExprMatrix[1:5,1:4]

# 2. SCENIC Analysis ----

# Load single-cell dataset
library(data.table)
melanomaData=fread("/Volumes/omics4tb/alomana/projects/mscni/data/testing.txt",sep="\t")

dataLength=dim(melanomaData)[1]
dataWidth=dim(melanomaData)[2]

geneNames=names(unlist(melanomaData[1,4:dataWidth]))
cellNames=unname(unlist(melanomaData[1:dataLength,1]))
melanomaExpression=t(as.matrix(melanomaData[1:dataLength,4:dataWidth]))
TPMplusOne=10**melanomaExpression
numericalError=10**melanomaExpression[1,1]
TPMs=TPMplusOne-numericalError

rownames(TPMs)=geneNames
colnames(TPMs)=cellNames

# calling SCENIC ----
dir.create("SCENIC_yapeng")
setwd("SCENIC_yapeng")

library(SCENIC)
org="hgnc" # or hgnc, or dmel
dbDir="/Volumes/omics4tb/alomana/projects/mscni/src/scenic/cisTarget_databases" # RcisTarget databases location
myDatasetTitle="SCENIC yapeng data analysis" # choose a name for your analysis
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=4) 

dir.create("int")

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Databases:
scenicOptions@settings$dbs <- setNames("hg19-tss-centered-10kb-7species.mc9nr.feather", "hg19-500bp-upstream-7species.mc9nr.feather")
scenicOptions@settings$db_mcVersion <- "v9"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# Analyze matrix for noise removal
nCellsPerGene <- apply(TPMs, 1, function(x) sum(x>0))
nCountsPerGene <- apply(TPMs, 1, sum)

summary(nCellsPerGene)
summary(nCountsPerGene)
max(TPMs)
sum(TPMs>0) / sum(TPMs==0)

# Filter 1: Remove cells that on average have fewer than 5 TPMs in 1% of cells
minReads <- 5*.01*ncol(TPMs)
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)]
length(genesLeft_minReads)

#Filter 2: Remove genes that appear in fewer than 1% of cells
minSamples <- ncol(TPMs)*.01
nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
length(genesLeft_minCells)

#Filter 3: Only keep the genes in the reference databases
library(RcisTarget)
motifRankings <- importRankings(getDatabases(scenicOptions)[[1]]) # either one, they should have the same genes
genesInDatabase <- colnames(getRanking(motifRankings))

working line

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

