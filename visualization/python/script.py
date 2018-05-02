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

def entropyCalculator(v,ground,sky):

    '''
    This function computes the entropy of a vector under defined range of expression
    '''

    # f.1. calculate the probability distribution
    step=(sky-ground)/len(cellIDs)
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
            clusterID=int(vector[1])
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

    geneNames.sort()
    cellIDs.sort()

    return expression,metadata,geneNames,cellIDs

def tSNERunner(task):

    '''
    This function calls t-SNE and makes a figure.
    '''

    print('\t\t working with task {}...'.format(task))
    
    thePerplexity=task[0]
    theLearningRate=task[1]

    # method='exact'

    embedded=sklearn.manifold.TSNE(perplexity=thePerplexity,learning_rate=theLearningRate,n_components=2,n_iter=10000,n_iter_without_progress=1000,init='pca',verbose=True).fit_transform(log10TPMsPO)

    figureName='figures.log10TPMsPO.barnes_hut/figure.tsne.p{}.lr{}.pdf'.format(thePerplexity,theLearningRate)
    matplotlib.pyplot.scatter(embedded[:,0],embedded[:,1],c=orderedColors,alpha=0.5,edgecolors='none')
    for i in range(len(selectedColors)):
        matplotlib.pyplot.scatter([],[],c=selectedColors[i],alpha=0.5,label=groupLabels[i],edgecolors='none')
    matplotlib.pyplot.legend()
    matplotlib.pyplot.xlabel('tSNE1')
    matplotlib.pyplot.ylabel('tSNE2')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()
    print()

    return None

###
### MAIN
###

# 0. user defined variables
dataFilePath='/Volumes/omics4tb/alomana/projects/malko/data/single.cell.data.for.Serdar.and.Adrian.txt'
numberOfThreads=4
selectedColors=['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
groupLabels=['State 1','State 2','State 3','State 4']

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

# correct expression from negative values and provide TPMs, log10 TPM+1 and log2 TPM+1
expressionMatrix=numpy.array(expressionMatrix)
TPMsPlusOne=10**expressionMatrix
TPMs=TPMsPlusOne-TPMsPlusOne[0,0]

log2TPMsPO=numpy.log2(TPMs+1)
log10TPMsPO=numpy.log10(TPMs+1)

ground=numpy.min(TPMs); sky=numpy.max(TPMs)
ground2=numpy.min(log2TPMsPO); sky2=numpy.max(log2TPMsPO)
ground10=numpy.min(log10TPMsPO); sky10=numpy.max(log10TPMsPO)

print('\t TPMs: ground {}; sky {}.'.format(ground,sky))
print('\t log2 (TPMs+1): ground {}; sky {}.'.format(ground2,sky2))
print('\t log10 (TPMs+1): ground {}; sky {}.'.format(ground10,sky10))

print('processing metadata...')
orderedAnnotation=[metadata[cellID] for cellID in cellIDs]
orderedColors=[selectedColors[annotation-1] for annotation in orderedAnnotation]

# 3. analyse data
print('analyzing data...')

"""

# 3.1. run a histogram of expression
print('\t building a histogram of expression...')
positiveValues=[]; nonPositiveValues=[]
for cell in log2TPMsPO:
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

# 3.2. build histogram of entropy per gene
print('\t building a histogram of entropy...')
entropies=[]
for gene in numpy.transpose(log2TPMsPO):
    s=entropyCalculator(gene,ground2,sky2)
    entropies.append(s)

x,y=histogrammer(entropies)

matplotlib.pyplot.plot(x,y,'-',color='black')

# compute entropy of references
ymax=[1/len(cellIDs) for element in range(len(cellIDs))]
smax=scipy.stats.entropy(ymax,base=2)
matplotlib.pyplot.axvline(smax,color='red',ls='--')
matplotlib.pyplot.axvline(0,color='blue',ls='--')

matplotlib.pyplot.xlabel('Gene entropy (bit)')
matplotlib.pyplot.ylabel('Probability')
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figure.entropy.distribution.pdf')
matplotlib.pyplot.clf()


# 3.3. run PCA and tSNE
print('\t dimensionality reduction...')

# 3.3.1. PCA
pca=sklearn.decomposition.PCA(n_components=2,svd_solver='full')
pca.fit(log2TPMsPO)
rotated=pca.transform(log2TPMsPO)

matplotlib.pyplot.scatter(rotated[:,0],rotated[:,1],c=orderedColors,alpha=0.5,edgecolors='none')

for i in range(len(selectedColors)):
    matplotlib.pyplot.scatter([],[],c=selectedColors[i],alpha=0.5,label=groupLabels[i],edgecolors='none')

matplotlib.pyplot.legend()
matplotlib.pyplot.xlabel('PC1')
matplotlib.pyplot.ylabel('PC2')
matplotlib.pyplot.savefig('figures/figure.pca.pdf')
matplotlib.pyplot.clf()

"""

# 3.3.2. t-SNE
perplexities=numpy.arange(5,50+5,5) # 5 to 30 seems the best range.
learningRates=numpy.arange(100,1000+100,100)

# consider putting it all into a single figure. use pickle to create all the embeddings. then do the plots
for thePerplexity in perplexities:
    for theLearningRate in learningRates:
        task=[thePerplexity,theLearningRate]
        tSNERunner(task)


# need to run with method='exact'

# consider making a simple heatmap with clustering for all genes, then only for the ones DETs between clusters. (remove zeros)

# define which genes are different between the different clusters, 5 by 5. Then represent a heatmap maybe?

### dist of max, min, mean, averages, number of zeros (droppouts)

    # consider making a simple heatmap with clustering
