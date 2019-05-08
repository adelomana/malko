import numpy,pandas,datetime,sklearn,s_dbw,sys
import matplotlib,matplotlib.pyplot
import multiprocessing,multiprocessing.pool
import scanpy
scanpy.settings.verbosity=5

def robustFinder(nei):

    '''
    This function returns the mean and standard deviation of goodness of fit for bootstrapped data.
    '''

    allK=[];allSS=[];allCHI=[];allVI=[]
    for iteration in range(bootstrapIterations):

        print(iteration)
        
        new=adata.copy()
        
        newIndices=numpy.random.choice(new.shape[0],new.shape[0])
        new.X=new.X[newIndices,:]
        
        scanpy.tl.pca(new,svd_solver='arpack')
        scanpy.pp.neighbors(new,n_neighbors=nei,n_pcs=50)
        scanpy.tl.umap(new)
        scanpy.tl.louvain(new)
        scanpy.pl.umap(new,color=['louvain'],palette='Set3',save='new.louvain.pdf',show=False)

        scanpy.tl.pca(adata,svd_solver='arpack')
        scanpy.pp.neighbors(adata,n_neighbors=nei,n_pcs=50)
        scanpy.tl.umap(adata)
        scanpy.tl.louvain(adata)
        scanpy.pl.umap(adata,color=['louvain'],palette='Set3',save='old.louvain.pdf',show=False)

        positions=new.obsm['X_umap']
        categories=new.obs['louvain'].tolist()
        numericCategories=[float(element) for element in categories]
        uniqueCategories=list(set(numericCategories))

        # goodness of partition
        SS=sklearn.metrics.silhouette_score(positions,numericCategories,metric='euclidean')
        CHI=sklearn.metrics.calinski_harabaz_score(positions,numericCategories)
        VI=s_dbw.S_Dbw(positions,numericCategories)

        print(len(uniqueCategories),SS,CHI,VI)
        
        allK.append(len(uniqueCategories)); allSS.append(SS); allCHI.append(CHI); allVI.append(VI)

    meanResult=(nei,numpy.mean(allK),numpy.mean(allSS),numpy.mean(allCHI),numpy.mean(allVI))
    stdResult=(nei,numpy.std(allK),numpy.std(allSS),numpy.std(allCHI),numpy.std(allVI))

    print(meanResult)

    return [meanResult,stdResult]

###
### MAIN
###

# 0. User defined variables
bootstrapIterations=3
numberOfNeighbors=numpy.arange(10,40+10,10)
print(numberOfNeighbors)
numberOfThreads=4

# 1. Read data
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t READ DATA FILE"))
#idata=scanpy.read_csv('/Volumes/omics4tb2/alomana/projects/mscni/data/scanpy/count.file.all.day.clean.csv')
idata=scanpy.read_csv('data/count.file.all.day.clean.csv')
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
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t FIND NUMBER OF CELL STATES"))
scanpy.settings.verbosity=0

#hydra=multiprocessing.pool.Pool(numberOfThreads)
#results=hydra.map(robustFinder,numberOfNeighbors)

results=[]
for nei in numberOfNeighbors:
    result=robustFinder()
    results.append(result)








# 3.3. save results, plot elsewhere.
#pickel!
    
# make sure that changing one single cell hardly changes anything. Otherwise, think about why?

print('\t generate figure...')


print(datetime.datetime.now().strftime("\t %Y-%m-%d %H:%M:%S"))
print('... completed.')
