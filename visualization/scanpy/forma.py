dataFolder='/Volumes/omics4tb2/alomana/projects/mscni/results/scanpy/resultsFile.12104/'

# read the gene names
dataFile=dataFolder+'var.csv'
geneNames=[]

with open(dataFile,'r') as f:
    next(f)
    for line in f:
        v=line.split(',')
        geneName=v[0]
        geneNames.append(geneName)

print(geneNames[:10],len(geneNames))

# read cell names
dataFile=resultsDir+'var.csv'
