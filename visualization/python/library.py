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
        for numberOfClusters in range(3,8+1):
            
            km=sklearn.cluster.KMeans(n_clusters=numberOfClusters, random_state=1).fit(embedded)
            kmLabels=km.labels_
            
            gmmLabels=sklearn.mixture.GaussianMixture(n_components=numberOfClusters,covariance_type='full').fit(embedded).predict(embedded)

            # f.3. compute goodness of clustering
            kmSS=sklearn.metrics.silhouette_score(embedded,kmLabels,metric='euclidean')
            gmmSS=sklearn.metrics.silhouette_score(embedded,gmmLabels,metric='euclidean')

            particularRanks.append(numberOfClusters); particularMethods.append('km'); particularQualifications.append(kmSS)
            particularRanks.append(numberOfClusters); particularMethods.append('gmm'); particularQualifications.append(gmmSS)


        # f.3. selecting best particular partition
        particularBestQualification=max(particularQualifications)
        particularBestRank=particularRanks[particularQualifications.index(particularBestQualification)]
        particularBestMethod=particularMethods[particularQualifications.index(particularBestQualification)]

        overallRanks.append(particularBestRank)
        overallMethods.append(particularBestMethod)
        overallQualifications.append(particularBestQualification)

    # f.4. selecting best overall partition
    a=max(overallQualifications)
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
