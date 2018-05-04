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

def generalAnalyzer(task):

    '''
    This function directs the main analysis.
    '''

    GOFs=[]
    for iteration in tsneRuns:
        
        # run  tSNEne
        thePerplexity=task[0]
        theLearningRate=task[1]
        tsneRunner(thePerplexity,theLearningRate)
        # perform clustering
    
        #compute goodness of clustering
    
    
    thePerplexity=task[0]
    theLearningRate=task[1]

    

    

    return None

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


def tsneRunner(thePerplexity,theLearningRate):

    '''
    This function runs tSNE.
    '''

    # method='exact'
    
    embedded=sklearn.manifold.TSNE(perplexity=thePerplexity,learning_rate=theLearningRate,n_components=2,n_iter=10000,n_iter_without_progress=1000,init='pca',verbose=True).fit_transform(TPMs)

    #figureName='figures.TPMsPO.barnes_hut/figure.tsne.p{}.lr{}.pdf'.format(thePerplexity,theLearningRate)
    #matplotlib.pyplot.scatter(embedded[:,0],embedded[:,1],c=orderedColors,alpha=0.5,edgecolors='none')
    #for i in range(len(selectedColors)):
    #    matplotlib.pyplot.scatter([],[],c=selectedColors[i],alpha=0.5,label=groupLabels[i],edgecolors='none')
    #matplotlib.pyplot.legend()
    #matplotlib.pyplot.xlabel('tSNE1')
    #matplotlib.pyplot.ylabel('tSNE2')
    #matplotlib.pyplot.tight_layout()
    #matplotlib.pyplot.savefig(figureName)
    #matplotlib.pyplot.clf()
    #print()

    return embedded

###
### MAIN
###

# 0. user defined variables
dataFilePath='/Volumes/omics4tb/alomana/projects/mscni/data/single.cell.data.txt'
numberOfThreads=4
selectedColors=['tab:blue', 'tab:green', 'tab:red', 'tab:purple']
groupLabels=['State 1','State 2','State 3','State 4']

#perplexities=numpy.arange(10,35+5,5) 
#learningRates=numpy.arange(100,800+100,100)
perplexities=numpy.arange(10,20+5,5) 
learningRates=numpy.arange(100,300+100,100)

tsneRuns=1 # this could be 3

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

# 3.0. define the parameters to test


tasks=[]; results=[]
for thePerplexity in perplexities:
    for theLearningRate in learningRates:
        task=[thePerplexity,theLearningRate]
        tasks.append(task)

print('\t {} tasks defined.'.format(len(tasks)))

# 3.1. sequential calling
for task in tasks:
    result=generalAnalyzer(taks)
    results.append(result)
    
# 3.2. parallel calling

# 4. plot the results
print(results)

# need to run with method='exact'

# consider making a simple heatmap with clustering for all genes, then only for the ones DETs between clusters. (remove zeros)

# define which genes are different between the different clusters, 5 by 5. Then represent a heatmap maybe?

### dist of max, min, mean, averages, number of zeros (droppouts)

    # consider making a simple heatmap with clustering
