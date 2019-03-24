###
### This script runs SCENIC for the complete trajectory of MSCNI.
###

# 0. SCENIC Set-up ----

# # 0.1. SCENIC dependencies ---
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GENIE3", "AUCell", "RcisTarget"), version = "3.8")
# # Also required:
# install.packages('zoo')

# # Recommended to run AUCell:
# BiocManager::install(c("mixtools", "rbokeh"))
# # To visualize the binary matrices and perform t-SNEs:
# BiocManager::install(c("NMF", "pheatmap", "Rtsne", "R2HTML"))
# # To support paralell execution (not available in Windows):
# BiocManager::install(c("doMC", "doRNG"))
# # To export/visualize in http://scope.aertslab.org
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

# # Other dependencies for the examples (lower priority)
# BiocManager::install(c("SingleCellExperiment"))

# # 0.2. install SCENIC ---
# # install.packages("devtools")
# devtools::install_github("aertslab/SCENIC", ref="v1.1.0")
# packageVersion("SCENIC")

# # 0.3. download human database, specify working directory first
# setwd("/Volumes/omics4tb2/alomana/projects/mscni/results/SCENIC/")
# dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
#              "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# # mc9nr: Motif collection version 9: 24k motifs
# 
# # dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
# for(featherURL in dbFiles)
# {
#   download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
#   descrURL <- gsub(".feather$", ".descr", featherURL)
#   if(file.exists(descrURL)) download.file(descrURL, destfile=basename(descrURL))
# }

# End SCENIC Set-up

# loading libraries
#### AUCell and RcisTarget were updated recently (May 2018/Jan 2019). Please, make sure you have the recommended versions: AUCell >=1.4.1 (minimum 1.2.4), RcisTarget>=1.2.0 (minimum 1.0.2), and GENIE3>=1.4.0 (minimum 1.2.1).
packageVersion("AUCell")
packageVersion("RcisTarget")
packageVersion("GENIE3") 

library('data.table')
library("GENIE3")
library("AUCell")
library("RcisTarget")
library('zoo')
library('tictoc')
library("mixtools")
library("rbokeh")
library("NMF")
library("pheatmap")
library("Rtsne")
library("R2HTML")
library("doMC")
library("doRNG")
library('devtools')
library('SingleCellExperiment')
library('SCENIC')

# 1. Set directories----
dir.create("/Volumes/omics4tb2/alomana/projects/mscni/results/SCENIC/results.2019.03.19")
setwd("/Volumes/omics4tb2/alomana/projects/mscni/results/SCENIC/results.2019.03.19") 

# 2. Read data----
melanomaData=fread("/Volumes/omics4tb2/alomana/projects/mscni/data/scanpy/count.file.all.day.clean.csv",sep=",")
dataLength=dim(melanomaData)[1]
dataWidth=dim(melanomaData)[2]
cellNames=names(unlist(melanomaData[1,2:dataWidth]))
length(cellNames)
head(cellNames)
geneNames=unname(unlist(melanomaData[1:dataLength,1]))
length(geneNames)
head(geneNames)
melanomaExpression=as.matrix(melanomaData[1:dataLength,2:dataWidth])
rownames(melanomaExpression)=geneNames
colnames(melanomaExpression)=cellNames
dim(melanomaExpression)

# 3. initialize scenic----
org="hgnc" # or hgnc, or dmel
dbDir="/Volumes/omics4tb2/alomana/projects/mscni/results/SCENIC/" # RcisTarget databases location
myDatasetTitle="SCENIC MSCNI 3k cells" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=8) 

# 4. filter data----
# (Adjust minimum values according to your dataset)
genesKept <- geneFiltering(melanomaExpression, scenicOptions=scenicOptions,minCountsPerGene=3*.01*ncol(melanomaExpression),minSamples=ncol(melanomaExpression)*.005)

interestingGenes <- c("SOX9", "SOX10", "DLX5")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- melanomaExpression[genesKept, ]
dim(exprMat_filtered)

# 5. correlation----
tic()
runCorrelation(exprMat_filtered, scenicOptions)
toc()

# 6. GENIE3----
# Optional: add log (if it is not logged/normalized already)
log2exprMat_filtered <- log2(exprMat_filtered+1) 
# Run GENIE3
tic()
runGenie3(log2exprMat_filtered, scenicOptions)
toc()

save.image(file ='/Volumes/omics4tb2/alomana/projects/mscni/results/SCENIC/tempo.7completed.2019.03.20.RData')

# 7. Build and score the GRN ----
load(file ='/Volumes/omics4tb2/alomana/projects/mscni/results/SCENIC/tempo.7completed.2019.03.20.RData',verbose = TRUE)
setwd("/Volumes/omics4tb2/alomana/projects/mscni/results/SCENIC/results.2019.03.19") 
scenicOptions@settings$nCores <- 2

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions) 
runSCENIC_3_scoreCells(scenicOptions, log2exprMat_filtered)

# 8. binarize network activity
runSCENIC_4_aucell_binarize(scenicOptions)

# 9. Clustering / dimensionality reduction on the regulon activity
nPcs <- c(5,15,50)
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5)
