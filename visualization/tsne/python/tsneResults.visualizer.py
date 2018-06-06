###
### This script visualizes the results from tsne.explorer.py.
###

import sys,numpy,pickle,seaborn,pandas
import sklearn,sklearn.cluster
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42


def embeddedGrapher():

    '''
    This function creates a figure for the embedded results.
    '''

    print('\t plotting figure {}/{}...'.format(i+1,len(sortedFitness)))

    embedded=embeddingBox[sortedFitness[i]][0]
    figurePerplexity=embeddingBox[sortedFitness[i]][1]
    figureLearningRate=embeddingBox[sortedFitness[i]][2]
    figureNC=embeddingBox[sortedFitness[i]][3]

    # f.1. clustering
    km=sklearn.cluster.KMeans(n_clusters=figureNC,random_state=1,algorithm='auto').fit(embedded)
    kmLabels=km.labels_
    kmCHI=numpy.log10(sklearn.metrics.calinski_harabaz_score(embedded,kmLabels))

    figureGF=kmCHI
    figureLabels=kmLabels
                
    # f.2. plot figure
    figureName='figures/figure.visual.{}.pdf'.format(str(i).zfill(5))
    cm=matplotlib.cm.get_cmap('tab20')
    colors50=list(cm.colors)+list(cm.colors)+list(cm.colors)[:10]
    orderedColors=[colors50[label-1] for label in figureLabels]

    matplotlib.pyplot.scatter(embedded[:,0],embedded[:,1],c=orderedColors,alpha=0.5,edgecolors='none')

    matplotlib.pyplot.title('CHI = {:0.3f}; nc = {}; P = {}; LR = {}'.format(figureGF,figureNC,figurePerplexity,figureLearningRate))
    matplotlib.pyplot.xlabel('tSNE1')
    matplotlib.pyplot.ylabel('tSNE2')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

    return None

# 0. user defined variables
similarityJar='/Volumes/omics4tb/alomana/projects/mscni/results/results.chi.iter100.2018.06.04.850.pickle'
figureName='/Volumes/omics4tb/alomana/projects/mscni/results/figures/figure.chi100.850.png'

scaleLimits=[3,3.4]
scaleIncrement=0.1
theMetricLabel='log$_{10}$ Calinski-Harabaz index'
minNC=3; maxNC=50

# 1. recover data
print('recovering data...')
f=open(similarityJar,'rb')
results=pickle.load(f)
f.close()

# 2. build heatmap figure
print('building heatmapfigure...')

# 2.1. recover parameter explored
SCs={};embeddingBox={}; ncDist={}

for result in results:
    perplexity=result[0][0]
    learningRate=result[0][1]
    SC=result[1]
    nc=result[2]
    method=result[3]
    print(method)
    embedding=result[4]

    # fill SC
    if perplexity not in SCs:
        SCs[perplexity]={}
    SCs[perplexity][learningRate]=[SC,nc]

    # fill embeddingBox
    embeddingBox[SC]=[embedding,perplexity,learningRate,nc]

    # fill ncDist
    if nc not in ncDist:
        ncDist[nc]=[SC]
    else:
        ncDist[nc].append(SC)

# 2.2. build a matrix of SCs
perplexities=list(SCs.keys())
learningRates=list(SCs[perplexities[0]].keys())
perplexities.sort(); learningRates.sort()

Mlist=[]; Nlist=[]
topLoc=[]; bestResult=0
for perplexity in perplexities:
    m=[]; n=[]
    for learningRate in learningRates:
        m.append(SCs[perplexity][learningRate][0]); n.append(SCs[perplexity][learningRate][1])
        if SCs[perplexity][learningRate][0] > bestResult:
            topLoc=[perplexity,learningRate]; bestResult=SCs[perplexity][learningRate][0]
    Mlist.append(m); Nlist.append(n)
M=numpy.array(Mlist); N=numpy.array(Nlist)
M=numpy.log10(M)

#print(M)
print('\t max value {}'.format(numpy.max(M)))
print('\t min value {}'.format(numpy.min(M)))
print('\t limits {}'.format(scaleLimits))
print('\t best results parameters {}'.format(topLoc))

# 2.3. build the figure
figureTitle='Goodness of clustering after tSNE'
matplotlib.pyplot.imshow(M,interpolation='none',cmap='viridis',vmin=scaleLimits[0],vmax=scaleLimits[1])
selectedTicks=list(numpy.arange(scaleLimits[0],scaleLimits[1]+scaleIncrement,scaleIncrement))
cb=matplotlib.pyplot.colorbar(orientation='vertical',fraction=0.05,ticks=selectedTicks)
cb.set_label(label=theMetricLabel,size=16)
cb.ax.tick_params(labelsize=16)
matplotlib.pyplot.grid(False)

# adding nc
x=-0.2
y=0.2
deltax=1.
deltay=1.
for i in range(len(N)):
    for j in range(len(N[0])):
        value=str(N[i][j])
        matplotlib.pyplot.text(x+deltax*j,y+deltay*i,value,fontsize=8,color='black')

xtickpositions=[]; xticknames=[]
ytickpositions=[]; yticknames=[]
for i in range(len(learningRates)):
    if i%5 == 0:
        xtickpositions.append(i); xticknames.append(learningRates[i])
for i in range(len(perplexities)):
    if i%5 == 0:
        ytickpositions.append(i); yticknames.append(perplexities[i])

matplotlib.pyplot.xticks(xtickpositions,xticknames,size=12,rotation=90)
matplotlib.pyplot.yticks(ytickpositions,yticknames,size=12)

matplotlib.pyplot.xlabel('Learning rate',fontsize=14)
matplotlib.pyplot.ylabel('Perplexity',fontsize=14)
       
matplotlib.pyplot.tick_params(axis='x',which='both',bottom='off',top='off')
matplotlib.pyplot.tick_params(axis='y',which='both',right='off',left='off')
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.title(figureTitle,fontsize=18)
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()

# 3. running embeddings
print('building embedding figures...')
sortedFitness=list(embeddingBox.keys())
sortedFitness.sort(reverse=True)

#for i in range(len(sortedFitness[:25])):
#    embeddedGrapher()

# 4. plot average fitness per nc

# 4.1. create a dataframe for plotting with seaborn
boxPlotPositions=[]
fitnessValues=[]

numberOfClusters=range(minNC,maxNC+1)
for nc in numberOfClusters:
    boxPlotPositions.append(nc)
    fitnessValues.append(None)

for nc in ncDist.keys():
    for value in ncDist[nc]:
        boxPlotPositions.append(nc)
        fitnessValues.append(numpy.log10(value))
# adding extra values
boxPlotPositions.append(5)
fitnessValues.append(None)

fitnessViolin=list(zip(boxPlotPositions,fitnessValues))
dfViolin=pandas.DataFrame(data=fitnessViolin,columns=['NC','Fitness'])

# 4.2 make figure
figureFile='figure.fitness.distribution.pdf'
#ax=seaborn.violinplot(x='NC',y='Fitness',data=dfViolin)
ax=seaborn.stripplot(x='NC',y='Fitness',data=dfViolin,jitter=True)

matplotlib.pyplot.ylim([3,3.4])
shownPositions=[3,5,7,10,13,20,30,40,50]
locations=list(numpy.array(shownPositions)-3)
matplotlib.pyplot.xticks(locations,shownPositions)
matplotlib.pyplot.grid(True,alpha=0.5,ls=':')
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(figureFile)
matplotlib.pyplot.clf()

# 4. last message
print('... analysis completed.')
