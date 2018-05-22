import numpy
import sklearn,sklearn.manifold,sklearn.cluster,sklearn.metrics,sklearn.mixture
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

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
dataFilePath='/Volumes/omics4tb/alomana/projects/mscni/data/single.cell.data.txt'

#selectedColors=['tab:blue', 'tab:green', 'tab:red', 'tab:purple', 'tab:orange', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
#groupLabels=['State 1','State 2','State 3','State 4']
thePerplexity=16
theLearningRate=950
trials=100
minNC=3; maxNC=20
numberOfClusters=numberOfClusters=range(minNC,maxNC+1)
figureName='figures/figure.tsne.p{}.lr{}.runner.chi.pdf'.format(thePerplexity,theLearningRate)

# 1. reading data
print('reading data...')
expression,metadata,geneNames,cellIDs=dataReader()
print('\t found {} cells with {} transcripts each.'.format(len(cellIDs),len(geneNames)))

print('processing metadata...')
orderedAnnotation=[metadata[cellID] for cellID in cellIDs]
cm=matplotlib.cm.get_cmap('tab20')
orderedColors=[cm.colors[annotation-1] for annotation in orderedAnnotation]

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
globalGF=0; globalLabels=[]

for trial in range(trials):
    print('\t working on trial {}'.format(trial))
    
    embedded=tsneRunner(thePerplexity,theLearningRate)

    localNC=0; localGF=0; locaLabels=None

    for nc in numberOfClusters:
            
            km=sklearn.cluster.KMeans(n_clusters=nc,random_state=1).fit(embedded)
            kmLabels=km.labels_
            
            # f.3. compute goodness of clustering
            #kmSS=sklearn.metrics.silhouette_score(embedded,kmLabels,metric='euclidean')
            kmCHI=numpy.log10(sklearn.metrics.calinski_harabaz_score(embedded,kmLabels))

            if kmCHI > localGF:
                localGF=kmCHI
                localNC=nc
                localLabels=kmLabels
                print('\t\t nc {} is best; GF {}'.format(nc,localGF))
                
    if localGF > globalGF and localGF < 3.8:
        globalGF=localGF
        finalLabels=localLabels

        print('\t\t\t improved GF {}'.format(globalGF))

        # 3.2. plott figure
        orderedColors=[cm.colors[label-1] for label in finalLabels]
        matplotlib.pyplot.scatter(embedded[:,0],embedded[:,1],c=orderedColors,alpha=0.5,edgecolors='none')

        matplotlib.pyplot.title('CHI = {}'.format(globalGF))

        matplotlib.pyplot.xlabel('tSNE1')
        matplotlib.pyplot.ylabel('tSNE2')
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig(figureName)
        matplotlib.pyplot.clf()

print('... analysis done.')

