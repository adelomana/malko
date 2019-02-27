###
### This script takes count data provided by collaborators and converts it into FPKM data for network inference.
###

# 0. incorporate Seurat
#install.packages('Seurat')
library(Seurat)

# 1. read data
expression=read.csv("/Volumes/omics4tb2/alomana/projects/mscni/data/seurat/count.file.alldays.csv",sep = ",",header=TRUE,row.names=1)

# 2. create a Seurat object and run a preliminary filter
seuratObject=CreateSeuratObject(raw.data=expression,min.cells=3,min.genes=200)

# 3. some plots about quality of run and filter out cells with too many UMIs
VlnPlot(object=seuratObject,features.plot=c("nGene","nUMI"),nCol=2)
GenePlot(object=seuratObject, gene1 = "nUMI", gene2 = "nGene")
seuratObject=FilterCells(object=seuratObject,subset.names=c("nUMI"),high.thresholds=c(4e4))

# 4. normalization
seuratObject=NormalizeData(object=seuratObject,normalization.method="LogNormalize",scale.factor=10000)

# 5. detection of varible genes
seuratObject=FindVariableGenes(object= seuratObject,mean.function=ExpMean,dispersion.function=LogVMR,x.low.cutoff=0.0125,x.high.cutoff=3,y.cutoff=0.5)
length(x=seuratObject@var.genes)

# 6. regress out unwanted sources of variation
seuratObject=ScaleData(object=seuratObject,vars.to.regress=c("nUMI"))

# head(seuratObject@raw.data)
# head(seuratObject@data)
# head(seuratObject@scale.data)

# 7. PCA
seuratObject=RunPCA(object=seuratObject,pc.genes=seuratObject@var.genes,do.print=TRUE,pcs.print=1:5,genes.print=5,pcs.compute=100)
VizPCA(object = seuratObject, pcs.use = 1:3)
PCAPlot(object = seuratObject, dim.1 = 1, dim.2 = 2)
seuratObject <- ProjectPCA(object = seuratObject, do.print = FALSE)
PCHeatmap(object = seuratObject, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = seuratObject, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
seuratObject=JackStraw(object=seuratObject,num.replicate=1000,display.progress=TRUE,num.pc=100,num.cores=8,do.par=TRUE)
JackStrawPlot(object=seuratObject,PCs=1:100)

PCElbowPlot(object=seuratObject,num.pc=50)

# 8. find clusters
seuratObject=FindClusters(object=seuratObject,reduction.type="pca",dims.use=1:40,resolution=0.6,print.output=0,save.SNN=TRUE)
PrintFindClustersParams(object=seuratObject)
table(seuratObject@ident)

# 9. tSNE visualization
seuratObject=RunTSNE(object=seuratObject,dims.use=1:40) 
TSNEPlot(object=seuratObject)

# 10. marker genes
markers1 = FindMarkers(seuratObject,1)
VlnPlot(object = seuratObject, features.plot = rownames(markers1)[1:5])
FeaturePlot(seuratObject, head(rownames(markers1)), cols.use = c("lightgrey", "blue"), nCol = 3)

markers <- FindAllMarkers(object = seuratObject, min.pct = 0.25, thresh.use = 0.25)
library(dplyr)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = seuratObject, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

# 10. export clean data
data_to_write_out <- as.data.frame(as.matrix(seuratObject@data))
library(data.table)
fwrite(x = data_to_write_out, row.names = TRUE, file = "/Volumes/omics4tb2/alomana/projects/mscni/results/seurat.cleaned.data.csv")

data_to_write_out <- as.data.frame(as.matrix(seuratObject@scale.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "/Volumes/omics4tb2/alomana/projects/mscni/results/seurat.cleaned.scaled.data.csv")

# 11. save work
saveRDS(seuratObject, file = "/Volumes/omics4tb2/alomana/projects/mscni/results/seurat.final.results.rds")