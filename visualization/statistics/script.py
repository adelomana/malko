###
### This script provides with the general distribution of expression and the number of dropouts
###

import sys,numpy

import matplotlib,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':24,'font.family':'Arial','xtick.labelsize':18,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def histogrammer(theData):

    '''
    This function creates a histogram.
    '''    

    x=[]; y=[]
    
    n,bins=numpy.histogram(theData,bins=int(numpy.sqrt(len(theData))))

    halfBin=(bins[1]-bins[0])/2.
    for bin in bins:
        center=bin+halfBin
        x.append(center)
    x.pop()
    y=numpy.array(n)
    y=list(y/float(sum(y)))

    return x,y

###
### MAIN
###

# 0. user defined variables
dataFilePath='/Volumes/omics4tb/alomana/projects/mscni/data/final/subset.txt'
dataFilePath='/Volumes/omics4tb/alomana/projects/mscni/data/final/subset2.txt'
dataFilePath='/Volumes/omics4tb/alomana/projects/mscni/data/final/Seurat_Scaled_Expression_M397_clean.txt'

# 1. read the file and provide with expression distribution and drop out rate
expression=[]
numberOfTranscripts=0

with open(dataFilePath,'r') as f:

    # obtain the number of cells
    header=f.readline()
    elements=header.split('\t')
    cellIDs=elements[1:]
    cellIDs[-1]=cellIDs[-1].replace('\n','')
    cellNumber=len(cellIDs)

    # obtain expression values
    for line in f:
        numberOfTranscripts=numberOfTranscripts+1
        v=line.split('\t')
        v[-1]=v[-1].replace('\n','')
        expressionValues=[float(element) for element in v[1:]]

        for element in expressionValues:
            expression.append(element)
            
# 2. make figures about expression distribution
print('number of cells: {}'.format(cellNumber))
print('number of transcripts: {}'.format(numberOfTranscripts))
print('max expression value: {}'.format(numpy.max(expression)))
print('min expression value: {}'.format(numpy.min(expression)))

# 3.1. run a histogram of expression
print('\t building a histogram of expression...')

sensibleExpression=[]
for element in expression:
    if element > 5:
        element = 5
    if element < -5:
        element = -5
    sensibleExpression.append(element)

x,y=histogrammer(sensibleExpression)

matplotlib.pyplot.plot(x,y,'-',color='black')
matplotlib.pyplot.grid(alpha=0.2)
matplotlib.pyplot.xlabel('log$_2$ scaled relative expression')
matplotlib.pyplot.ylabel('Probability')
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figure.expression.distribution.png',dpi=300)
matplotlib.pyplot.clf()
