pathToDirectory <- "/Users/MattWall/Desktop/network_package/data/Day.3.network_inference"
setwd(pathToDirectory)

# 1. Install SCENIC Dependencies ------------------------------------------

source("https://bioconductor.org/biocLite.R")
biocLite(c("GENIE3", "RcisTarget", "AUCell","data.table","SingleCellExperiment","tidyverse", "factoextra"))

#source("https://bioconductor.org/biocLite.R")
#biocLite(c("feather"))

# Recommended to run AUCell:
biocLite(c("mixtools"))
# To visualize the binary matrices and perform t-SNEs:
biocLite(c("NMF", "Rtsne", "R2HTML"))
# To support paralell execution:
biocLite(c("doMC", "doRNG"))

#install.packages("doMC", repos="http://R-Forge.R-project.org") #for pc
#library(doMC) #for pc

# compatability with notebooks
install.packages("knitr")

## Install SCENIC
install.packages("devtools")
devtools::install_github("aertslab/SCENIC")
install.packages('zoo')
install.packages("rbokeh")
install.packages("gplots")
install.packages("Cairo")
install.packages("stringr")

# pseudoinverse for transcription factor activity
install.packages("corpcor")

install.packages("readr")
library(readr)

# library with fast manipulations for dataframes
library(data.table)
# library to create a bioconductor object for SCENIC
library(SingleCellExperiment)
# library for motif enrichment (transcription factor binding site analysis)
library(RcisTarget)
# library with SCENIC functions for network analysis
library(SCENIC)
# libraries for t-SNE plotting
library(gplots)
library(KernSmooth)
library(RColorBrewer)
#library with pseudoinverse for transcription factor activity
library(corpcor)
library(stringr)


# Check versions. Correct versions are: AUCell 1.2.4, RcisTarget 1.0.2, and GENIE3 1.2.1 or posterior.
packageVersion("AUCell")
packageVersion("RcisTarget")
packageVersion("GENIE3")

##If you have trouble installing these packages with biocLite, try direct installation from website:
#install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/AUCell_1.2.4.tar.gz", repos=NULL)
#install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/RcisTarget_1.0.2.tar.gz", repos=NULL)
#install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/GENIE3_1.2.1.tar.gz", repos=NULL)

# Suppress loading messages when building the HTML
suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(SingleCellExperiment)
})

# Define functions for enrichment analysis

read_gmt <- function(filename){
  myList <- list()
  myFile <- read_lines(file=filename)
  for(i in 1:length(myFile)){
    tmpLine <- str_split(myFile[[i]],"\t")
    myList[[tmpLine[[1]][[1]]]] <- tmpLine[[1]][3:length(tmpLine[[1]])]
  }
  return(myList)
}

invert.list <- function(myList){
  invList <- list()
  for(element in names(myList)){
    tmpLine <- myList[[element]]
    for(g in (1:length(tmpLine))){
      gene <- tmpLine[[g]]
      if(gene %in% names(invList)){
        invList[[gene]] <- c(invList[[gene]],element)
      } else {
        invList[[gene]] <- c(element)
      }
    }
  }
  return(invList)
}

write.list <- function(myList,filename){
  sink(filename)
  nms <- names(myList)
  for(i in (1:(length(myList)-1))){
    line <- c(nms[[i]],myList[[i]])
    writeLines(line,sep = " ")
    writeLines("\n",sep = "")
  }
  line <- c(nms[[length(myList)]],myList[[length(myList)]])
  writeLines(line,sep = " ")
  sink()
}

read.list <- function(filename){
  myList <- list()
  myFile <- read_lines(file=filename)
  for(i in 1:length(myFile)){
    tmpLine <- str_split(myFile[[i]]," ")
    myList[[tmpLine[[1]][[1]]]] <- tmpLine[[1]][2:(length(tmpLine[[1]])-1)]
  }
  return(myList)
}

Enrichment <- function(Regulons,db,inv.db){
  
  db.genes <- names(inv.db)
  regulonNames <- names(Regulons)
  resultsList <- list()  
  for(regName in regulonNames){
    regGenes = Regulons[[regName]]
    common.genes <- intersect(regGenes,db.genes)
    Hits <- vector()
    for(name in common.genes){
      Hits <- c(Hits,inv.db[[name]])
    }
    pHtable <- table(Hits)
    testers <- pHtable[pHtable>2]
    if(length(testers) < 1) next
    pop <- length(db.genes)
    pList <- list()
    for(i in 1:length(testers)){
      m <- length(db[[names(testers)[[i]]]])
      draw <- length(common.genes)
      overlap <- testers[[i]]
      p <- phyper(overlap,m,pop-m,draw,lower.tail = F)
      pList[names(testers)[[i]]] <- p
    }
    resultsList[[regName]] <- pList
  }
  return(resultsList)
}

ttestByPatient <- function(patientLabel){
  patientCells <- rownames(cellLabels)[cellLabels[,1]==patientLabel]
  patientActivities <- tfActivities[,patientCells]
  othersActivities <- tfActivities[,setdiff(colnames(tfActivities),patientCells)]
  orderedCols <- c(patientCells,setdiff(colnames(tfActivities),patientCells))
  
  tfs <- rownames(patientActivities)
  ts <- vector()
  ps <- vector()
  for(i in 1:length(tfs)){
    ttest <- t.test(patientActivities[i,],othersActivities[i,])
    ts <- c(ts,ttest$statistic)
    ps <- c(ps,ttest$p.value)
  }
  sortTtest <- order(ts,decreasing = T)
  tt.results <- list()
  tt.results[["regulon"]] <- tfs[sortTtest]
  tt.results[["t statistic"]] <- ts[sortTtest]
  tt.results[["p-value"]] <- ps[sortTtest]
  tt.df <- as.data.frame(tt.results)
  return(tt.df)
}

# 1.1. Import databases -----------------------------------------------------

hallmarksFile <- "functional_databases/hallmarks.v6.2.symbols.gmt"
pathwaysFile <- "functional_databases/pathways.v6.2.symbols.gmt"
goBPFile <- "functional_databases/goBP.v6.2.symbols.gmt"
goCCFile <- "functional_databases/goCC.v6.2.symbols.gmt"
goMFFile <- "functional_databases/goMF.v6.2.symbols.gmt"

# read in the database files
hallmarks <- read_gmt(hallmarksFile)
pathways <- read_gmt(pathwaysFile)
goBP <- read_gmt(goBPFile)
goCC <- read_gmt(goCCFile)
goMF <- read_gmt(goMFFile)

inverseHallmarks <- invert.list(hallmarks)
inversePathways <- invert.list(pathways)
inverseGOBP <- invert.list(goBP)
inverseGOCC <- invert.list(goCC)
inverseGOMF <- invert.list(goMF)

## Download motif databases for RcisTarget (slow: can take >30 min)

# mc9nr: Motif collection version 9: 24k motifs
#dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather","https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")

#dir.create("RcisTarget")
#setwd("RcisTarget") # if needed
#for(featherURL in dbFiles)
#{
#  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
#  descrURL <- gsub(".feather$", ".descr", featherURL)
#  if(file.exists(descrURL)) download.file(descrURL, destfile=basename(descrURL))
#}

# End SCENIC Set-up -------------------------------------------------------



# SCENIC Analysis ---------------------------------------------------------

# 2. Create Directory to Save SCENIC Files --------------------------------

setwd(pathToDirectory)
#dir.create("SCENIC_Melanoma") # for future use, create folder to save files
setwd("SCENIC_Melanoma") # set to folder created in previous line

# 3. Load Files for SCENIC Run --------------------------------------------

# 3.1 Load expression file
expressionFile <- "../input_data/malignant.8kgenes.data.csv" #"../" is how we code the "back" arrow to go up one folder
#expressionFile <- "../input_data/mm.8kgenes.data.csv"

scData <- data.table::fread(expressionFile, sep=",")
geneNames <- unname(unlist(scData[,1, with=FALSE]))
exprMatrix <- as.matrix(scData[,-1, with=FALSE])
rm(scData)
dim(exprMatrix)
rownames(exprMatrix) <- geneNames

# Check that matrix is labeled properly and shows expression values
exprMatrix[1:5,1:4]

# Remove duplicated rows
exprMatrix <- exprMatrix[unique(rownames(exprMatrix)),] 
dim(exprMatrix) # number of genes and number of samples

# 3.2 Load cell Labels file
# Label cell types
labelsFile <- "../input_data/malignant.8kgenes.tumorMetadata.csv" # "../" is how we write the "back" arrow to go up one folder

scMetadata <- read.csv(labelsFile, row.names=1,header=TRUE)
colnames(scMetadata) <- "CellType"
cellLabels <- as.data.frame(scMetadata)

# 3.3 Create single Bioconductor object using SingleCellExperiment

scHumanMelanoma <- SingleCellExperiment(assays = list(counts = exprMatrix),
                                      colData=data.frame(cellLabels[colnames(exprMatrix),, drop=FALSE]))

# save scHumanMelanoma bioconductor object
save(scHumanMelanoma, file="SCENIC_data/scHumanMelanoma.RData")

# End Data Loading for SCENIC Run -----------------------------------------


# 4. Checkpoint to Start SCENIC from Files --------------------------------

# 4.1 Load bioconductor object
load("SCENIC_data/scHumanMelanoma.RData")
exprMat <- counts(scHumanMelanoma)
dim(exprMat)

# 4.2 Collect information used in plot labels
cellInfo <- colData(scHumanMelanoma)
cellInfo$nGene <- colSums(exprMat>0)
cellInfo <- data.frame(cellInfo)
head(cellInfo)
#dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

# 4.3 update colVars for new dataset
# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=setNames(c("forestgreen", "darkorange", "magenta4", "hotpink", "red3", "skyblue"), 
                                  c("Mel78", "Mel79", "Mel80", "Mel81", "Mel88", "Mel89")))
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

# 4.4 Set inputs for SCENIC run

org="hgnc" # or hgnc, or dmel
dbDir="RcisTarget" # RcisTarget databases location
myDatasetTitle="SCENIC example on Human Melanoma" # choose a name for your analysis
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=4) 

# Modify for plotting if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Add database location to scenicOptions:
scenicOptions@settings$dbs <- setNames("hg19-tss-centered-10kb-7species.mc9nr.feather", "hg19-500bp-upstream-7species.mc9nr.feather")
scenicOptions@settings$db_mcVersion <- "v9"

# Save to scenicOptions to use at a later time
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


# 5. Filter Dataset for Analysis ------------------------------------------

# 5.1 Analyze matrix for noise removal
nCellsPerGene <- apply(exprMat, 1, function(x) sum(x>0))
nCountsPerGene <- apply(exprMat, 1, sum)

summary(nCellsPerGene)
summary(nCountsPerGene)
max(exprMat)
sum(exprMat>0) / sum(exprMat==0)

# Filter 1: Remove cells that on average have fewer than 1 log2(TPM+1) in 1% of cells

minReads <- 1*.01*ncol(exprMat)
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)]
length(genesLeft_minReads)

#Filter 2: Remove genes that appear in fewer than 1% of cells
minSamples <- ncol(exprMat)*.01
nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
length(genesLeft_minCells)

#Filter 3: Only keep the genes in the reference databases

#motifRankings <- importRankings(getDatabases(scenicOptions)[[1]]) # either one, they should have the same genes
#genesInDatabase <- colnames(getRanking(motifRankings))
#genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)] # Uncomment when doing future analyses

# read file for today's session. Use commented function above for future sessions
#melanoma file
melanoma_genesLeft_minCells_inDatabases <- read.csv("genesLeft_minCells_inDatabases.csv") # This loads the results from line 300. Don't use this line in future work.
genesLeft_minCells_inDatabases <- as.character(melanoma_genesLeft_minCells_inDatabases["x"][[1]])
length(genesLeft_minCells_inDatabases)


# Check whether any relevant gene / potential gene of interest is missing:
#interestingGenes <- c("AHR","ZFX")
#interestingGenes[which(!interestingGenes %in% genesLeft_minCells_inDatabases)]

genesKept <- genesLeft_minCells_inDatabases
#saveRDS(genesKept, file=getIntName(scenicOptions, "genesKept"))

# This is your filtered expression matrix
exprMat_filtered <- exprMat[genesKept, ]

# To avoid confusion in the coming steps:
rm(exprMat)

# 6. Run SCENIC -----------------------------------------------------------

# 6.1 Generate correlation matrix
## Perform Correlation
#corrMat <- cor(t(exprMat_filtered), method="spearman")
## (Only the rows for TFs will be needed):
#allTFs <- getDbTfs(scenicOptions)
#corrMat <- corrMat[which(rownames(corrMat) %in% allTFs),]
#saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))

# 6.2 Run GENIE3 (Do not do this today!)

#runGenie3(exprMat_filtered, scenicOptions)

# 6.3 Load bioconductor object into Matrix form for SCENIC analysis

load("SCENIC_data/scHumanMelanoma.RData")
exprMat <- counts(scHumanMelanoma)
dim(exprMat)

#Set SCENIC parameters
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 4
scenicOptions@settings$seed <- 123

# 6.4 Run SCENIC using wrappers (Do not do this today)
#runSCENIC_1_coexNetwork2modules(scenicOptions)
#runSCENIC_2_createRegulons(scenicOptions)
#runSCENIC_3_scoreCells(scenicOptions, exprMat)

# End SCENIC run ----------------------------------------------------------


# 7. Visualize results with t-SNE -----------------------------------------

#Set SCENIC parameters

# 7.1 Grid search of different t-SNE settings (takes about 2 minutes):
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(20,30,40), perpl=c(20,30,40))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(20,30,40), perpl=c(20,30,40), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot gridsearch t-SNE results as pdf (individual files in int/): 
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE)

# 7.2 Plot t-SNE standard results:
par(mfcol=c(3,3))
par(mar=c(1,1,1,1)) #This prevents the plotting error "Error in plot.new() : figure margins too large"
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

# 7.3 Using only "high-confidence" regulons (normally similar)
par(mfcol=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

# 7.4 Set tSNE defaults 
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 30
scenicOptions@settings$defaultTsne$perpl <- 40
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

## Option to binarize regulons as "on" or "off"
# scenicOptions@settings$devType="png"
#runSCENIC_4_aucell_binarize(scenicOptions)

# 7.5 Optimize plots with App

logMat <- exprMat # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat) #default t-SNE
savedSelections <- shiny::runApp(aucellApp)


# 7.6 Generate composite plot of all regulons:

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

tfsToPlot <- c("XBP1", "ETV1", "TP53")
#regulonPdfFile <- "output/regulonPlots"
#pdf(paste(paste(regulonPdfFile,"expression",sep = "_"),"pdf",sep = "."))
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[tfsToPlot],], plots="Expression")
#dev.off()

# Save all AUC into one PDF:
Cairo::CairoPDF("output/Composite_tSNE_plots.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()


# 8. Explore Regulons -----------------------------------------------------

regulons <- loadInt(scenicOptions, "regulons")

# 8.1 open a list of regulons for inspection
regulons[c("XBP1","ETV1")]

#open a single regulon
regulonName <- "XBP1"
regulonGenes <- regulons[c(regulonName)][[1]]

# 8.2 create incidence matrix - important for later use
incidList <- reshape2::melt(regulons)
incidMat <- table(incidList[,2], incidList[,1])
incidMat[1:5,1:5]
saveRDS(incidMat, file=getIntName(scenicOptions, "regulons_incidMat"))

splitnames <- str_split_fixed(row.names(incidMat),"_", n=2)
conciseIncidMat <- incidMat[splitnames[,2] !="extended",]
extendedIncidMat <- incidMat[splitnames[,2] =="extended",]
conciseIncidMat[1:5,1:5]

# 8.3 Search and display regulon motifs
#regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
#regulonTargetsInfo <- RcisTarget::addLogo(regulonTargetsInfo, motifCol="bestMotif")
#regulonTargetsInfo$Genie3Weight <- signif(regulonTargetsInfo$Genie3Weight, 2)

#install.packages("DT")
#colsToShow <- c("TF", "gene", "nMotifs", "bestMotif", "logo", "NES", "highConfAnnot", "Genie3Weight")
#DT::datatable(regulonTargetsInfo[TF=="ETV1" & highConfAnnot==TRUE, colsToShow, with=F], escape=FALSE, filter="top")


# 9. Compute Regulator Activities -----------------------------------------

tfCoefficients <- t(conciseIncidMat)
reducedExprMat <- exprMat_filtered[colnames(conciseIncidMat),]
concisePseudo <- pseudoinverse(tfCoefficients)
tfActivities <- concisePseudo %*% reducedExprMat
row.names(tfActivities) <- row.names(conciseIncidMat)

#perform ttest of regulon activities of your patient vs. the others
patient <- "Mel80"
tt.df <- ttestByPatient(patientLabel=patient)

#plot first 5 rows of ttest dataframe
tt.df[1:5,]

#save the ttest dataframe as a csv file
write.csv(tt.df,paste(paste("output/ttest",patient,sep = "_"),"csv",sep = "."),row.names=FALSE)

#plot the activities of the First regulon

regulonLabel <- "XBP1"
patientCells <- rownames(cellLabels)[cellLabels[,1]==patient]
patientActivities <- tfActivities[,patientCells]
othersActivities <- tfActivities[,setdiff(colnames(tfActivities),patientCells)]
orderedCols <- c(patientCells,setdiff(colnames(tfActivities),patientCells))
plot(tfActivities[regulonLabel,orderedCols],col=cellLabels[orderedCols,1],xlab="Patient index",ylab="Regulon activity")

#save plot to pdf file
regulonActivityPdfFile <- paste("output",paste(patient,regulonLabel,sep = "_"),sep="/")
pdf(paste(paste(regulonActivityPdfFile,"activity",sep = "_"),"pdf",sep = "."))
plot(tfActivities[regulonLabel,orderedCols],col=cellLabels[orderedCols,1],xlab="Patient index",ylab="Regulon activity")
dev.off()

# 10. Functional analysis -------------------------------------------------

# 10.1. Compute Enrichments -----------------------------------------------------

hallmarkEnrichment <- Enrichment(regulons,hallmarks,inverseHallmarks)
pathwayEnrichment <- Enrichment(regulons,pathways,inversePathways)
GOBPEnrichment <- Enrichment(regulons,goBP,inverseGOBP)
GOCCEnrichment <- Enrichment(regulons,goCC,inverseGOCC)
GOMFEnrichment <- Enrichment(regulons,goMF,inverseGOMF)

# 10.2. View Results ---------------------------------------------------------

# set resultsList to the functional enrichment you wish to review
resultsList <- GOBPEnrichment
regulonHits <- names(resultsList)

# List of regulons that have at least one possible hit
#print(regulonHits)

# Check if regulon had functional enrichment for selected database (output = TRUE/FALSE)

regulonName <- "XBP1"
regulonName %in% regulonHits

# Generate dataframe of results for regulon
hitResults <- resultsList[regulonName]
viewResults <- t(as.data.frame(hitResults)[order(as.data.frame(hitResults))])
colnames(viewResults) <- "p"

# View only the top hits
#numberHits <- 50
#numberToView <- min(numberHits,dim(viewResults)[[1]]) 
#head(viewResults,numberToView)

# View only hits with BH-corrected p-value < 0.05
correctedP <- p.adjust(viewResults[,"p"], method = "BH", n = length(viewResults))
correctedP <- as.data.frame(correctedP)
colnames(correctedP) <- "corrected p-value"
if(dim(correctedP)[1] == 1){
  viewResults
} else head(correctedP,sum(correctedP[,"corrected p-value"]<0.05))

# write the results of the functional analysis

write.csv(correctedP,paste(paste(paste("output/Enrichment",patient,sep = "_"),regulonName,sep="_"),"csv",sep="."))

# look at the genes in a pathway
functionalTerm <- "REACTOME_UNFOLDED_PROTEIN_RESPONSE"
functionalGenes <- pathways[[functionalTerm]]
functionalGenes
# look at the genes in the regulon 

regulonGenes # defined above in section 8.1

# look at the genes in the regulon that appear in the functional gene set
intersect(functionalGenes,regulonGenes)
