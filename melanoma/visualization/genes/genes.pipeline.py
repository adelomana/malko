import numpy,pandas,datetime,sklearn,s_dbw,sys,pickle
import multiprocessing,multiprocessing.pool
import scanpy
scanpy.settings.verbosity=5

def robustFinder(nei):

    '''
    This function returns the mean and standard deviation of goodness of fit for bootstrapped data.
    '''
    
    # work on original data
    scanpy.tl.pca(adata,svd_solver='arpack')
    scanpy.pp.neighbors(adata,n_neighbors=nei,n_pcs=50)
    scanpy.tl.umap(adata)
    scanpy.tl.louvain(adata)

    positions=adata.obsm['X_pca']
    categories=adata.obs['louvain'].tolist()
    numericCategories=[float(element) for element in categories]
    uniqueCategories=list(set(numericCategories))

    # goodness of partition
    SS=sklearn.metrics.silhouette_score(positions,numericCategories,metric='euclidean')
    CHI=sklearn.metrics.calinski_harabaz_score(positions,numericCategories)
    VI=s_dbw.S_Dbw(positions,numericCategories)
    original=(nei,len(uniqueCategories),SS,CHI,VI)
    
    allK=[];allSS=[];allCHI=[];allVI=[]
    for iteration in range(bootstrapIterations):

        # changes
        new=adata.copy()
        allIndices=numpy.arange(adata.shape[0])
        deemed=numpy.random.randint(adata.shape[0])
        allIndicesButOne=numpy.delete(allIndices,deemed)
        new=new[allIndicesButOne,:]
        
        scanpy.tl.pca(new,svd_solver='arpack')
        scanpy.pp.neighbors(new,n_neighbors=nei,n_pcs=50)
        scanpy.tl.umap(new)
        scanpy.tl.louvain(new)
        #scanpy.pl.umap(new,color=['louvain'],palette='Set3',save='new.nei{}.iter{}.louvain.pdf'.format(nei,iteration),show=False)
        #matplotlib.pyplot.close()

        positions=new.obsm['X_pca']
        categories=new.obs['louvain'].tolist()
        numericCategories=[float(element) for element in categories]
        uniqueCategories=list(set(numericCategories))

        # goodness of partition
        SS=sklearn.metrics.silhouette_score(positions,numericCategories,metric='euclidean')
        CHI=sklearn.metrics.calinski_harabaz_score(positions,numericCategories)
        VI=s_dbw.S_Dbw(positions,numericCategories)

        #print('neighbors: {}; iter: {}; clusters: {}; SS: {:.2E}; CHI: {:.2E}; VI: {:.2E}'.format(nei,iteration,len(uniqueCategories),SS,CHI,VI))
        
        allK.append(len(uniqueCategories)); allSS.append(SS); allCHI.append(CHI); allVI.append(VI)

    meanResult=(nei,numpy.mean(allK),numpy.mean(allSS),numpy.mean(allCHI),numpy.mean(allVI))
    stdResult=(nei,numpy.std(allK),numpy.std(allSS),numpy.std(allCHI),numpy.std(allVI))

    return [original,meanResult,stdResult]

###
### MAIN
###

# 0. User defined variables
bootstrapIterations=25
numberOfNeighbors=numpy.arange(5,40+5,5)
print(numberOfNeighbors)
numberOfThreads=2
jar='jars/exploration.results.003.pickle'

# 1. Read data
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t READ DATA FILE"))
idata=scanpy.read_csv('/Volumes/omics4tb2/alomana/projects/mscni/data/scanpy/count.file.all.day.clean.csv')
#idata=scanpy.read_csv('data/count.file.all.day.clean.csv')
adata=idata.transpose()
print(adata)

# 2. Preprocessing
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t PREPROCESSING"))

scanpy.pp.filter_cells(adata,min_genes=200)
scanpy.pp.filter_genes(adata,min_cells=3)

adata.obs['n_counts'] = adata.X.sum(axis=1)

scanpy.pp.normalize_per_cell(adata, counts_per_cell_after=1e5)
scanpy.pp.log1p(adata)

adata.raw = adata

scanpy.pp.highly_variable_genes(adata,min_mean=0.0125,max_mean=6,min_disp=0.25) # 2,851
adata = adata[:, adata.var['highly_variable']]

scanpy.pp.regress_out(adata, ['n_counts'])
scanpy.pp.scale(adata, max_value=10)

# 3. Find number of states
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t EXPLORE PARAMETER SPACE"))
scanpy.settings.verbosity=0

#hydra=multiprocessing.pool.Pool(numberOfThreads)
#results=hydra.map(robustFinder,numberOfNeighbors)

results=[]
for nei in numberOfNeighbors:
    result=robustFinder(nei)
    results.append(result)

# 4. save results
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t SAVE RESULTS"))
f=open(jar,'wb')
pickle.dump(results,f)
f.close()
    
# 5. final call
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t COMPLETED"))
