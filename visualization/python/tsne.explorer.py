###
### This script visualizes in an optimal manner melanoma single-cell transcriptomes.
###

import sys,numpy,pickle
import scipy,scipy.stats
import sklearn,sklearn.manifold,sklearn.cluster,sklearn.metrics,sklearn.mixture
import multiprocessing,multiprocessing.pool

import library

###
### MAIN
###

# 0. user defined variables
dataFilePath='/Volumes/omics4tb/alomana/projects/mscni/data/single.cell.data.txt'
numberOfThreads=4

perplexities=numpy.arange(10,35+5,5) 
learningRates=numpy.arange(100,800+100,100)

tsneRuns=10 

# 1. reading data
print('reading data...')
expression,metadata,geneNames,cellIDs=library.dataReader()
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
#for task in tasks:
#    result=generalAnalyzer(task)
#    results.append(result)
    
# 3.1.b. parallel calling
hydra=multiprocessing.pool.Pool(numberOfThreads)
results=hydra.map(generalAnalyzer,tasks)

# 3.2. pickle the results
similarityJar='results.pickle'
f=open(similarityJar,'wb')
pickle.dump(results,f)
f.close()
