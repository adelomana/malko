import sys,numpy

dataFolder='/Volumes/omics4tb2/alomana/projects/mscni/results/scanpy/resultsFile.1273/'

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
dataFile=dataFolder+'obs.csv'
cellIDs=[]

with open(dataFile,'r') as f:
    next(f)
    for line in f:
        v=line.split(',')
        cellID=v[0]+'.'+v[-1].replace('\n','')
        cellIDs.append(cellID)

print(cellIDs[:10],len(cellIDs))

# read expression values
dataFile=dataFolder+'X.csv'
expression=[]
with open(dataFile,'r') as f:
    for line in f:
        v=line.split(',')
        values=[float(element) for element in v]
        
        expression.append(values)

print(len(expression))
for i in range(20):
    tpmsPerCell=[10**element for element in expression[i]]
    print(min(tpmsPerCell),max(tpmsPerCell),numpy.log10(sum(tpmsPerCell)))

# writing a file
formattedFile=dataFolder+'formatted.expression.data.csv'

with open(formattedFile,'w') as f:
    for cellID in cellIDs:
        f.write(',{}'.format(cellID))
    f.write('\n')
    
    for i in range(len(geneNames)):
        f.write('{}'.format(geneNames[i]))
        for j in range(len(cellIDs)):
            f.write(',{}'.format(expression[j][i]))
        f.write('\n')
