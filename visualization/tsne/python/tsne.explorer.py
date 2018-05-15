###
### This script visualizes in an optimal manner melanoma single-cell transcriptomes.
###

import sys,numpy,pickle,time
import scipy,scipy.stats
import sklearn,sklearn.manifold,sklearn.cluster,sklearn.metrics,sklearn.mixture
import multiprocessing,multiprocessing.pool

import matplotlib,matplotlib.pyplot

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

    print('\t working with task {}'.format(task))

    overallRanks=[]; overallMethods=[]; overallQualifications=[]
    
    for iteration in range(tsneRuns):
        
        # f.1. run  tSNEne
        thePerplexity=task[0]
        theLearningRate=task[1]
        embedded=tsneRunner(thePerplexity,theLearningRate)

        # f.2. perform clustering and goodness of clustering
        particularRanks=[]; particularMethods=[]; particularQualifications=[]
        numberOfClusters=range(minNC,maxNC+1)
        for nc in numberOfClusters:
            
            km=sklearn.cluster.KMeans(n_clusters=nc,random_state=1).fit(embedded)
            kmLabels=km.labels_
            
            gmmLabels=sklearn.mixture.GaussianMixture(n_components=nc,covariance_type='full').fit(embedded).predict(embedded)

            # f.3. compute goodness of clustering
            kmSS=sklearn.metrics.silhouette_score(embedded,kmLabels,metric='euclidean')
            gmmSS=sklearn.metrics.silhouette_score(embedded,gmmLabels,metric='euclidean')

            kmCHI=sklearn.metrics.calinski_harabaz_score(embedded,kmLabels)
            gmCHI=sklearn.metrics.calinski_harabaz_score(embedded,gmmLabels)

            #particularRanks.append(nc); particularMethods.append('km'); particularQualifications.append(kmSS)
            #particularRanks.append(nc); particularMethods.append('gmm'); particularQualifications.append(gmmSS)

            particularRanks.append(nc); particularMethods.append('km'); particularQualifications.append(kmCHI)
            particularRanks.append(nc); particularMethods.append('gmm'); particularQualifications.append(gmCHI)

        #particularCHI=list(numpy.array(particularCHI)/max(particularCHI))
        #matplotlib.pyplot.plot(range(minNC,maxNC+1),particularSC,'ok',label='SC')
        #matplotlib.pyplot.plot(range(minNC,maxNC+1),particularCHI,'or',label='CHI')
        #matplotlib.pyplot.legend()
        #localName='figures/figure.{}.pdf'.format(time.time())
        #matplotlib.pyplot.savefig(localName)
        #matplotlib.pyplot.clf()

        #bestCHI=particularCHI.index(max(particularCHI))
        #bestSC=particularSC.index(max(particularSC))
                                        
        #print('best SC {}; best CHI {}'.format(numberOfClusters[bestSC],numberOfClusters[bestCHI]))

        # f.3. selecting best particular partition
        particularBestQualification=max(particularQualifications)
        particularBestRank=particularRanks[particularQualifications.index(particularBestQualification)]
        particularBestMethod=particularMethods[particularQualifications.index(particularBestQualification)]

        overallRanks.append(particularBestRank)
        overallMethods.append(particularBestMethod)
        overallQualifications.append(particularBestQualification)

    # f.4. selecting best overall partition
    a=max(overallQualifications)
    #a=numpy.median(overallQualifications)
    b=overallRanks[overallQualifications.index(a)]
    c=overallMethods[overallQualifications.index(a)]

    result=[task,a,b,c]
    
    return result

def tsneRunner(thePerplexity,theLearningRate):

    '''
    This function runs tSNE.
    '''

    # method='exact'
    
    embedded=sklearn.manifold.TSNE(perplexity=thePerplexity,learning_rate=theLearningRate,n_components=2,n_iter=10000,n_iter_without_progress=1000,init='pca',verbose=0).fit_transform(log2TPMsPO)
    
    return embedded

###
### MAIN
###

# 0. user defined variables
dataFilePath='/proj/omics4tb/alomana/projects/mscni/data/single.cell.data.txt'
resultsJar='results.10.50.100.800.bestnc.tsne25runs.pickle'
numberOfThreads=16
perplexities=numpy.arange(10,50+2,2) 
learningRates=numpy.arange(100,800+100,100)
tsneRuns=25
minNC=3; maxNC=25

# 1. reading data
print('reading data...')
expression,metadata,geneNames,cellIDs=dataReader()
print('\t found {} cells with {} transcripts each.'.format(len(cellIDs),len(geneNames)))

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

# 3. analyse data
print('analyzing data...')

# 3.0. define the parameters to test
tasks=[]; results=[]
for thePerplexity in perplexities:
    for theLearningRate in learningRates:
        task=[thePerplexity,theLearningRate]
        tasks.append(task)
print('\t {} tasks defined.'.format(len(tasks)))

# 3.1.a. sequential calling
for task in tasks:
    result=generalAnalyzer(task)
    results.append(result)

# 3.1.b. parallel calling
#hydra=multiprocessing.pool.Pool(numberOfThreads)
#results=hydra.map(generalAnalyzer,tasks)

# 3.2. pickle the results
f=open(resultsJar,'wb')
pickle.dump(results,f)
f.close()
