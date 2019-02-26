# 0. SCENIC Set-up ----

# 0.1. SCENIC dependencies ---
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
library("mixtools")
# To visualize the binary matrices and perform t-SNEs:
library("NMF")
library("Rtsne")
library("R2HTML")
# To support paralell execution:
library("doMC")
library("doRNG")

library("RCurl")
library("XML")
library(data.table)
library('zoo')

# To visualize in http://scope.aertslab.org
#install.packages("devtools")
library(devtools)
install_github("aertslab/SCopeLoomR")

# 0.2. Install SCENIC ----
setwd("/Volumes/omics4tb/alomana/projects/mscni/results/scenic.2019/")

# install.packages("devtools")
install_github("aertslab/SCENIC")
library(SCENIC)

# 0.3. Download databases ----
## Download motif databases for RcisTarget (slow: can take >30 min)

dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs

# dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
  descrURL <- gsub(".feather$", ".descr", featherURL)
  if(file.exists(descrURL)) download.file(descrURL, destfile=basename(descrURL))
}
# End SCENIC Set-up

# 1. Melanoma data loading and filtering ----

# 1.1. Load single-cell dataset ----
melanomaData=fread("/Volumes/omics4tb/alomana/projects/mscni/data/single.cell.data.txt",sep="\t")

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

# 1.2. filter data for SCENIC ----
dir.create("SCENIC_one.day_2019.01.23")
setwd("SCENIC_one.day_2019.01.23")

org="hgnc" # or hgnc, or dmel
dbDir="/Volumes/omics4tb/alomana/projects/mscni/results/scenic.2019/cisTarget_databases" # RcisTarget databases location
myDatasetTitle="SCENIC.2019.01.23.one.day" # choose a name for your analysis
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=8) 

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

genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)

genesKept <- genesLeft_minCells_inDatabases
saveRDS(genesKept, file=getIntName(scenicOptions, "genesKept"))

TPMs_filtered <- TPMs[genesKept, ]
rm(TPMs)

## 1.3. Perform Correlation ---
corrMat <- cor(t(TPMs_filtered), method="spearman")
# (Only the rows for TFs will be needed needed):
# allTFs <- getDbTfs(scenicOptions)
# corrMat <- corrMat[which(rownames(corrMat) %in% allTFs),]
saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))


# 2. Run SCENIC --------------------------------------------------------------

# Run GENIE3
# setwd("SCENIC_MouseBrain")
# library(SCENIC)
# scenicOptions <- readRDS("int/scenicOptions.Rds")
# library(SingleCellExperiment)
# load("data/sceMouseBrain.RData")
# exprMat <- counts(sceMouseBrain)
# genesKept <- loadInt(scenicOptions, "genesKept")
# exprMat_filtered <- exprMat[genesKept,]

# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(TPMs_filtered+1) 

# 2.1. Run GENIE3 ----
library(tictoc)
tic()
runGenie3(exprMat_filtered, scenicOptions)
toc()

# 2.2. Build and score the GRN (runSCENIC) ----
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 8
scenicOptions@settings$seed <- 123

tic()
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
toc()

save.image('results.up2.section3.line174.RData')

# 3. visualize results ----
scenicOptions@settings$seed <- 123 # same seed for all of them
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(5,25,50), perpl=c(10,15,20,25))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(5,25,50), perpl=c(10,15,20,25), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/): 
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE)

par(mfcol=c(3,4))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5)

# Using only "high-confidence" regulons (normally similar)
par(mfcol=c(3,4))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5)

scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 50
scenicOptions@settings$defaultTsne$perpl <- 25
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# binarization of regulons
library('rbokeh')
logMat <- exprMat_filtered # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# scenicOptions@settings$devType="png"
exprMat <- exprMat_filtered 
runSCENIC_4_aucell_binarize(scenicOptions)

# 4. exploring results ----
logMat <- exprMat # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat) #default t-SNE
savedSelections <- shiny::runApp(aucellApp)

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("ATF4", "HOXB2", "MITF","HES1", "RXRG", "BRCA1")],], plots="Expression")

save.image('up2line227.RData')

# 8. Explore Regulons -----------------------------------------------------

lst <- loadInt(scenicOptions, "regulons")

max_l <- max(rapply(lst, length))
df=as.data.frame(rapply(lst, function(x) 'length<-'(x, max_l), how="list"))
tdf=transpose(df)
colnames(tdf) <- rownames(df)
rownames(tdf) <- colnames(df)
write.csv(tdf,file='testing.csv')

