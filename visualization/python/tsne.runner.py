import numpy
import matplotlib,matplotlib.pyplot

import library

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

###
### MAIN
###

# 0. user defined variables
dataFilePath='/Volumes/omics4tb/alomana/projects/mscni/data/single.cell.data.txt'

selectedColors=['tab:blue', 'tab:green', 'tab:red', 'tab:purple', 'tab:orange', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
groupLabels=['State 1','State 2','State 3','State 4']

# 1. reading data
print('reading data...')
expression,metadata,geneNames,cellIDs=library.dataReader(dataFilePath)
print('\t found {} cells with {} transcripts each.'.format(len(cellIDs),len(geneNames)))

print('processing metadata...')
orderedAnnotation=[metadata[cellID] for cellID in cellIDs]
orderedColors=[selectedColors[annotation-1] for annotation in orderedAnnotation]

# 2. process data
print('processing data...')
expressionMatrix=[]
for cellID in cellIDs:
    row=[]
    for geneName in geneNames:
        value=expression[cellID][geneName]
        row.append(value)
    expressionMatrix.append(row)
expressionMatrix=numpy.array(expressionMatrix)

# correct expression from negative values and provide TPMs, log10 TPM+1 and log2 TPM+1
expressionMatrix=numpy.array(expressionMatrix)
TPMsPlusOne=10**expressionMatrix
TPMs=TPMsPlusOne-TPMsPlusOne[0,0]

log2TPMsPO=numpy.log2(TPMs+1)

ground2=numpy.min(log2TPMsPO); sky2=numpy.max(log2TPMsPO)

print('\t log2 (TPMs+1): ground {}; sky {}.'.format(ground2,sky2))

# 3. analyse data
print('analyzing data...')

# 3.0. run tSNE
thePerplexity=25
theLearningRate=100
embedded=library.tsneRunner(thePerplexity,theLearningRate,log2TPMsPO)

# 3.2. plott figure
figureName='figures/figure.tsne.p{}.lr{}.runner.pdf'.format(thePerplexity,theLearningRate)
orderedColors=[selectedColors[label] for label in kmLabels]

matplotlib.pyplot.scatter(embedded[:,0],embedded[:,1],c=orderedColors,alpha=0.5,edgecolors='none')

matplotlib.pyplot.xlabel('tSNE1')
matplotlib.pyplot.ylabel('tSNE2')
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()

matplotlib.pyplot.clf()

print('... analysis done.')

