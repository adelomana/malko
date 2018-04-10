###
### This script analyses melanoma single-cell transcriptomes.
###

import sys,numpy
import scipy,scipy.stats
import sklearn,sklearn.decomposition,sklearn.manifold
import multiprocessing,multiprocessing.pool
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def entropyCalculator(v):

    '''
    This function computes the entropy of a vector under defined range of expression
    '''

    # f.1. calculate the probability distribution
    step=0.1
    k=numpy.arange(ground,sky+step,step)
    n,bins=numpy.histogram(v,bins=k)
    y=[]
    y=numpy.array(n)
    y=y/float(sum(y))

    # f.2. calculate entropy
    s=scipy.stats.entropy(y,base=2)

    return s


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

def dataReader():

    '''
    This function reads the expression file.
    '''

    expression={}
    metadata={}
    cellIDs=[]
    clusterIDs=[]

    with open(dataFilePath, 'r') as f:
        # work with header file
        header=f.readline()
        h=header.split('\t')
        h[-1]=h[-1].replace('\n','')
        geneNames=h[2:]
        # work with body file
        for line in f:
            vector=line.split('\t')
            cellID=vector[0]
            cellIDs.append(cellID)
            clusterID=[1]
            if clusterID not in clusterIDs:
                clusterIDs.append(clusterID)
            e=vector[2:]
            E=[float(element) for element in e]

            if len(E) != len(geneNames):
                print('error at reading expression...')
                sys.exit()

            # filling up variables
            expression[cellID]={}
            for i in range(len(E)):
                if geneNames[i] not in expression[cellID]:
                    expression[cellID][geneNames[i]]=E[i]

            metadata[cellID]=clusterID

    # check for redundancy
    #print(len(geneNames),len(cellIDs))
    #geneNames=list(set(geneNames))
    #cellIDs=list(set(cellIDs))
    #print(len(geneNames),len(cellIDs))

    geneNames.sort()
    cellIDs.sort()

    return expression,metadata,geneNames,cellIDs

def tSNERunner(task):

    '''
    This function calls t-SNE and makes a figure.
    '''

    thePerplexity=task[0]
    theLearningRate=task[1]

    embedded=sklearn.manifold.TSNE(perplexity=thePerplexity,learning_rate=theLearningRate,n_components=2,n_iter=10000,n_iter_without_progress=1000,init='pca',method='barnes_hut').fit_transform(expressionMatrix)

    figureName='figures/figure.tsne.p{}.lr{}.pdf'.format(thePerplexity,theLearningRate)
    
    matplotlib.pyplot.scatter(embedded[:,0],embedded[:,1])
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

    return None

###
### MAIN
###

# 0. user defined variables
dataFilePath='/Volumes/omics4tb/alomana/projects/malko/data/single.cell.data.for.Serdar.and.Adrian.txt'
numberOfThreads=-1

# 1. reading data
print('reading data...')
expression,metadata,geneNames,cellIDs=dataReader()
print('found {} cells with {} transcripts each.'.format(len(cellIDs),len(geneNames)))

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

# correct expression to non-negative data
print('\t removing negative values...')
expressionMatrix[expressionMatrix<0]=0
ground=numpy.min(expressionMatrix)
sky=numpy.max(expressionMatrix)
print('\t ground {}; sky {}.'.format(ground,sky))

# 3. analyse data
print('analyzing data...')

"""

# 3.1. run a histogram of expression
print('\t building a histogram of expression...')
positiveValues=[]; nonPositiveValues=[]
for cell in expressionMatrix:
    for value in cell:
        if value > 0:
            positiveValues.append(value)
        else:
            nonPositiveValues.append(value)

a=len(positiveValues); b=len(nonPositiveValues); c=a+b

print('\t {} values detected, {} positive, {} non-positive.'.format(c,a/c,b/c))

x,y=histogrammer(positiveValues)

matplotlib.pyplot.plot(x,y,'-',color='black')
matplotlib.pyplot.xlabel('Relative expression')
matplotlib.pyplot.ylabel('Probability')
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figure.expression.distribution.pdf')
matplotlib.pyplot.clf()

# 3.2. build histogram of entropy per cell
print('\t building a histogram of entropy...')
entropies=[]
for cell in expressionMatrix:
    s=entropyCalculator(cell)
    entropies.append(s)

x,y=histogrammer(entropies)

matplotlib.pyplot.plot(x,y,'-',color='black')
matplotlib.pyplot.xlabel('Entropy (bit)')
matplotlib.pyplot.ylabel('Probability')
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figure.entropy.distribution.pdf')
matplotlib.pyplot.clf()

"""

# 3.3. run PCA and tSNE
print('\t dimensionality reduction...')

"""

# 3.3.1. PCA
pca=sklearn.decomposition.PCA(n_components=2)
pca.fit(expressionMatrix)
rotated=pca.transform(expressionMatrix)

matplotlib.pyplot.scatter(rotated[:,0],rotated[:,1])
matplotlib.pyplot.savefig('figures/figure.pca.pdf')
matplotlib.pyplot.clf()

"""

# 3.3.2. t-SNE
perplexities=numpy.arange(5,50+5,5)
learningRates=numpy.arange(10,1000+10,10)

tasks=[]

for thePerplexity in perplexities:
    for theLearningRate in learningRates:
        tasks.append([thePerplexity,theLearningRate])


# parallel version
#hydra=multiprocessing.pool.Pool(numberOfThreads)
#tempo=hydra.map(tSNERunner,tasks)

# single-thread

for task in tasks:
    tSNERunner(task)
    sys.exit()

# consider making a simple heatmap with clustering for all genes, then only for the ones DETs between clusters. (remove zeros)

# define which genes are different between the different clusters, 5 by 5. Then represent a heatmap maybe?

### dist of max, min, mean, averages, number of zeros (droppouts)

    # consider making a simple heatmap with clustering
