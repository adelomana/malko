import sys,scanpy,pickle,numpy,time,sklearn,s_dbw

def preliminaryAnalysis():

    # f1. Read data
    idata=scanpy.read_csv('/Volumes/omics4tb2/alomana/projects/mscni/data/scanpy/count.file.all.day.clean.csv')
    adata=idata.transpose()

    # f2. Preprocessing
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

    scanpy.tl.pca(adata,svd_solver='arpack') ### there seem to be a bug

    return adata

def stateFinder(adata):

    # exploration
    results=[]

    scanpy.tl.pca(adata,svd_solver='arpack')
    
    for nei in numberOfNeighbors:

        # work on original data
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

        result=(nei,len(uniqueCategories),SS,CHI,VI)
        results.append(result)

    # storage
    jar='jars/jar.iter{}.pickle'.format(iteration)
    f=open(jar,'wb')
    pickle.dump(results,f)
    f.close()
    
    return 

###
### MAIN
###

# 0. user defined variables
scanpy.settings.verbosity=0
fine=numpy.arange(5,31,1)
coarse=numpy.arange(35,55,5)
veryCoarse=numpy.arange(75,225,25)
joined=numpy.concatenate((fine,coarse))
numberOfNeighbors=numpy.concatenate((joined,veryCoarse))
#numberOfNeighbors=numpy.arange(5,40+5,5)

# 1. reading iteration 
if sys.argv[1] == '-i':
    iteration=int(sys.argv[2])
else:
    raise ValueError('wrong flag')

if iteration < 0:
    original=True
else:
    original=False

# 2. run preliminary analysis
adata=preliminaryAnalysis()

# 3. complete analysis
if original == True:
    stateFinder(adata)
else:
    allIndices=numpy.arange(adata.shape[0])
    deemed=numpy.random.randint(adata.shape[0])
    allIndicesButOne=numpy.delete(allIndices,deemed)
    adata=adata[allIndicesButOne,:]
    stateFinder(adata)
