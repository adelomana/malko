import numpy,pandas,datetime,sys
import matplotlib,matplotlib.pyplot
import scanpy
scanpy.settings.verbosity=5

# 1. read data
print('reading data...')
adata=scanpy.read_10x_mtx('/Volumes/omics4tb2/alomana/projects/gbm/results/with/outs/filtered_feature_bc_matrix',var_names='gene_symbols',cache=True)
adata.var_names_make_unique() 
print(adata)
print()

# 2. QC
print('before QC...')
print(adata)

print('QC...')
scanpy.pp.filter_cells(adata,min_genes=500)
scanpy.pp.filter_genes(adata,min_cells=10)
scanpy.pl.highest_expr_genes(adata,n_top=20,show=False,save='highestAbundant.pdf')

mitoGenes=adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = numpy.sum(adata[:, mitoGenes].X, axis=1).A1 / numpy.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

scanpy.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True,show=False,save='QC.violin.multipanel.pdf')
scanpy.pl.scatter(adata, x='n_counts', y='percent_mito',show=False,save='QC.counts.mito.0.pdf')
scanpy.pl.scatter(adata, x='n_counts', y='n_genes',show=False,save='QC.counts.genes.0.pdf')

print('remove high mito')
adata = adata[adata.obs['percent_mito'] < 0.05, :]
print(adata)
print()

scanpy.pl.scatter(adata, x='n_counts', y='percent_mito',show=False,save='QC.counts.mito.1.pdf')
scanpy.pl.scatter(adata, x='n_counts', y='n_genes',show=False,save='QC.counts.genes.1.pdf')

meanCount=numpy.mean(adata.obs['n_counts'])
standardDeviationCount=numpy.std(adata.obs['n_counts'])
suggestedRange=[meanCount-(2*standardDeviationCount),meanCount+(2*standardDeviationCount)]
print('mean count, standard deviation, suggested range')
print(meanCount,standardDeviationCount,suggestedRange)
print()

print('remove top count outliers')
adata = adata[adata.obs['n_counts'] < suggestedRange[1], :]
print(adata)
print()

print('remove bottom bottom outliers')
#adata = adata[adata.obs['n_counts'] > suggestedRange[0], :]
adata = adata[adata.obs['n_counts'] > 5000, :]
print(adata)
print()

print('remove low n_genes')
adata = adata[adata.obs['n_genes'] > 2000, :]
print(adata)
print()

scanpy.pl.scatter(adata, x='n_counts', y='percent_mito',show=False,save='QC.counts.mito.2.pdf')
scanpy.pl.scatter(adata, x='n_counts', y='n_genes',show=False,save='QC.counts.genes.2.pdf')

# 3. normalization
print('normalization')
scanpy.pp.normalize_per_cell(adata, counts_per_cell_after=suggestedRange[1])
scanpy.pp.log1p(adata)
adata.raw = adata

# 4. selection of highly-variable genes
#scanpy.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)  # n_obs × n_vars = 2673 × 1822
#scanpy.pp.highly_variable_genes(adata, min_mean=0.1, max_mean=4, min_disp=0.5)     # n_obs × n_vars = 2673 × 1386
#scanpy.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=4, min_disp=1/3)   # n_obs × n_vars = 2673 × 2430
scanpy.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=4, min_disp=0.5)   # n_obs × n_vars = 2673 × 1855 

scanpy.pl.highly_variable_genes(adata,show=False,save='highly.variable.genes.pdf')
print(adata)
adata = adata[:,adata.var['highly_variable']]
print(adata)

# 5. regress out and scale
scanpy.pp.regress_out(adata, ['n_counts', 'percent_mito']) # dont forget to regress bc of cell cycle
scanpy.pp.scale(adata, max_value=10)

# 6. associate sample metadata
print('associate sample metadata...')
cellIDs=adata.obs_names.tolist()

sampleLabels=[]
for cellID in cellIDs:
    if '-1' in cellID: 
        sampleLabels.append('T1')
    elif '-2' in cellID:
        sampleLabels.append('T2')
    elif '-3' in cellID:
        sampleLabels.append('T3')
    elif '-4' in cellID:
        sampleLabels.append('T4a')
    elif '-5' in cellID:
        sampleLabels.append('T4b')
    else:
        raise ValueError('cellID not recognized')
    
adata.obs['metadata']=sampleLabels
uniqueSampleLabels=list(set(sampleLabels))
uniqueSampleLabels.sort()
print()

# 7. PCA
print('PCA...')
scanpy.tl.pca(adata, svd_solver='arpack')
scanpy.pl.pca(adata, color='metadata',palette=['gray','blue','green','orange','red'],show=False,save='PCA.pdf')
scanpy.pl.pca(adata, color='FABP7',palette='viridis',show=False,save='PCA.FABP7.pdf')
scanpy.pl.pca_variance_ratio(adata, log=True,show=False,save='PCA.components.pdf',n_pcs = 50)
print()

# 8. tSNE
print('tSNE...')
scanpy.tl.tsne(adata)
scanpy.pl.tsne(adata,color='metadata',palette=['gray','blue','green','orange','red'],show=False,save='tSNE.pdf')
scanpy.pl.tsne(adata,color='n_counts',palette='viridis',show=False,save='tSNE.counts.pdf')
scanpy.pl.tsne(adata,color='n_genes',palette='viridis',show=False,save='tSNE.genes.pdf')
print()

# 9. UMAP resolution

"""
print('UMAP resolution...')
possibleNeighbors=numpy.arange(30,101,1)
louvainClusters=[]
print(possibleNeighbors)
for nei in possibleNeighbors:
    scanpy.pp.neighbors(adata,n_neighbors=nei,n_pcs=50,knn=True)

    scanpy.tl.umap(adata)
    scanpy.pl.umap(adata, color='metadata',palette=['blue','green','orange','red','gray'],show=False,save='umap.{}.pdf'.format(nei))
    scanpy.pl.umap(adata, color='n_counts',palette='viridis',show=False,save='umap.counts.{}.pdf'.format(nei))
    scanpy.pl.umap(adata, color='n_genes',palette='viridis',show=False,save='umap.genes.{}.pdf'.format(nei))

    # 10. Louvain
    scanpy.tl.louvain(adata)
    scanpy.pl.umap(adata, color='louvain',palette='tab10',show=False,save='umap.louvain.{}.pdf'.format(nei))

    uniqueClusters=adata.obs['louvain'].nunique()
    louvainClusters.append(uniqueClusters)

matplotlib.pyplot.plot(possibleNeighbors,louvainClusters,'o-',color='black')
matplotlib.pyplot.xlabel('neighbors')
matplotlib.pyplot.ylabel('clusters')
matplotlib.pyplot.savefig('figures/louvainRanks.pdf')
"""

# 10. UMAP & Louvain once 
print('UMAP single day...')
nei=90
scanpy.pp.neighbors(adata,n_neighbors=nei,n_pcs=50,knn=True)

scanpy.tl.umap(adata)
scanpy.pl.umap(adata, color='metadata',palette=['gray','blue','green','orange','red'],show=False,save='umap.{}.pdf'.format(nei))
scanpy.pl.umap(adata, color='n_counts',palette='viridis',show=False,save='umap.counts.{}.pdf'.format(nei))
scanpy.pl.umap(adata, color='n_genes',palette='viridis',show=False,save='umap.genes.{}.pdf'.format(nei))
scanpy.pl.umap(adata, color='FABP7',palette='viridis',show=False,save='umap.FABP7.{}.pdf'.format(nei))

scanpy.tl.louvain(adata)
scanpy.pl.umap(adata, color='louvain',palette='tab10',show=False,save='umap.louvain.{}.pdf'.format(nei))
print()

# 11. per day
print('UMAP per day...')
LouvainColors=adata.uns['louvain_colors']
LouvainMemberships=adata.obs['louvain'].tolist()
positions=adata.obsm['X_umap']
accumulatedDays=[]

for dayLabel in uniqueSampleLabels[:-1]:
    
    accumulatedDays.append(dayLabel)
    
    if dayLabel == 'T4a':
        accumulatedDays.append(uniqueSampleLabels[-1])

    tag='.'.join(accumulatedDays)
    
    xpos=[]; ypos=[]; myColors=[]
    for i in range(len(adata.obs['metadata'])):
        if adata.obs['metadata'][i] in accumulatedDays:
            xpos.append(positions[i,0])
            ypos.append(positions[i,1])
            myColors.append(LouvainColors[int(LouvainMemberships[i])])
    
    print(tag,accumulatedDays,len(xpos))
    matplotlib.pyplot.scatter(xpos,ypos,color=myColors,alpha=2/3,edgecolors='none')
    
    matplotlib.pyplot.xlim([-5,13])
    matplotlib.pyplot.ylim([-8.5,0.5])
    matplotlib.pyplot.xticks([])
    matplotlib.pyplot.yticks([])
    
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('figures.perDay/{}.perDay.louvain.pdf'.format(tag))
    matplotlib.pyplot.clf()
print('')

# 12. per day non-accumulative
print('non-accumulative...')
for dayLabel in uniqueSampleLabels[:-1]:
    
    accumulatedDays=[]
    accumulatedDays.append(dayLabel)
    if dayLabel == 'T4a':
        accumulatedDays.append(uniqueSampleLabels[-1])
        
    tag='.'.join(accumulatedDays)
    
    xpos=[]; ypos=[]; myColors=[]
    for i in range(len(adata.obs['metadata'])):
        if adata.obs['metadata'][i] in accumulatedDays:
            xpos.append(positions[i,0])
            ypos.append(positions[i,1])
            myColors.append(LouvainColors[int(LouvainMemberships[i])])
    
    print(tag,accumulatedDays,len(xpos))
    matplotlib.pyplot.scatter(xpos,ypos,color=myColors,alpha=2/3,edgecolors='none')
    
    matplotlib.pyplot.xlim([-5,13])
    matplotlib.pyplot.ylim([-8.5,0.5])
    matplotlib.pyplot.xticks([])
    matplotlib.pyplot.yticks([])
    
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('figures.perDay/{}.perDay.non.cumulative.louvain.pdf'.format(tag))
    matplotlib.pyplot.clf()
