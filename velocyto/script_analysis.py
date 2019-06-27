import sys,numpy,datetime,velocyto,matplotlib,logging

LOG_FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=LOG_FORMAT,level=logging.INFO,datefmt='%Y-%m-%d %H:%M:%S \t')

# 1. read data
logging.info('reading velocyto data')
vlm = velocyto.VelocytoLoom("/Volumes/omics4tb2/alomana/projects/mscni/results/velocyto/combined.loom")


# 2. read scanpy data


# 2. normalize
logging.info('normalization velocyto output')
vlm.normalize("S", size=True, log=True)
vlm.S_norm  # contains log normalized

#vlm.plot_fractions()
#matplotlib.pyplot.savefig('proportions.pdf')
#matplotlib.pyplot.clf()

# 3. preliminary filtering
vlm.filter_cells(bool_array=vlm.initial_Ucell_size > numpy.percentile(vlm.initial_Ucell_size, 0.5))

# 4. associate cluster id
logging.info('associate cluster id')
print(vlm,type(vlm))


sys.exit()
# 1. Read data
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t READ DATA FILE"))
idata=scanpy.read_csv('/Volumes/omics4tb2/alomana/projects/mscni/data/scanpy/count.file.all.day.clean.csv')
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
scanpy.tl.pca(adata,svd_solver='arpack')
scanpy.pp.neighbors(adata,n_neighbors=12,n_pcs=50)
scanpy.tl.umap(adata)
scanpy.tl.louvain(adata)
scanpy.pl.umap(adata,color=['louvain'],palette='Set3',save='new.nei.iter.louvain.pdf',show=False)
matplotlib.pyplot.clf()

##### read velocyto run
vlm = velocyto.VelocytoLoom("/Volumes/omics4tb2/alomana/projects/mscni/results/velocyto/18324_Yapeng_single_cell/velocyto/18324_Yapeng_single_cell.loom")

vlm.normalize("S", size=True, log=True)
vlm.S_norm  # contains log normalized


vlm.plot_fractions()
matplotlib.pyplot.savefig('proportions.pdf')
matplotlib.pyplot.clf()


vlm.filter_cells(bool_array=vlm.initial_Ucell_size > numpy.percentile(vlm.initial_Ucell_size, 0.5))


# vlm.set_clusters(vlm.ca["ClusterName"])

vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)
vlm.filter_genes(by_detection_levels=True)


vlm.score_cv_vs_mean(3000, plot='filter.pdf', max_expr_avg=35)
vlm.filter_genes(by_cv_vs_mean=True)

vlm._normalize_S(relative_size=vlm.S.sum(0),
             target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0),
             target_size=vlm.U.sum(0).mean())

### 5 prep


vlm.perform_PCA()
nei=int(0.025*969)
print(nei)
vlm.knn_imputation(k=nei,n_pca_dims=20,balanced=True,n_jobs=8)
vlm.fit_gammas()


vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)


from sklearn.manifold import TSNE
bh_tsne = TSNE()
vlm.ts = bh_tsne.fit_transform(vlm.pcs[:, :25])


vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,n_neighbors=nei, knn_random=True, sampled_fraction=0.5)
vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)
vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=nei)


matplotlib.pyplot.scatter(vlm.embedding[:, 0], vlm.embedding[:, 1],c="0.8", alpha=0.2, s=10, edgecolor="")
matplotlib.pyplot.savefig('figure.simple.pdf')
matplotlib.pyplot.clf()


# overlay arrows
ix_choice = numpy.random.choice(vlm.embedding.shape[0], size=int(vlm.embedding.shape[0]/1.), replace=False)
matplotlib.pyplot.scatter(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],c="0.8", alpha=0.4, s=10, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)

quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,linewidths=0.25, width=0.00045,edgecolors="k", color='black', alpha=1)
matplotlib.pyplot.quiver(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],vlm.delta_embedding[ix_choice, 0], vlm.delta_embedding[ix_choice, 1],**quiver_kwargs)

matplotlib.pyplot.axis("off")
matplotlib.pyplot.savefig("day.24.full_arrows.pdf")
matplotlib.pyplot.clf()
