###
### This script visualizes the results from tsne.explorer.py.
###

import sys,numpy,pickle
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

# 0. user defined variables
similarityJar='results.pickle'

# 1. recover data
print('recovering data...')
f=open(similarityJar,'rb')
results=pickle.load(f)
f.close()

# 2. build figure
print('building figure...')

# 2.1. recover parameter explored
SCs={}

for result in results:
    perplexity=result[0][0]
    learningRate=result[0][1]
    SC=result[1]
    nc=result[2]

    if perplexity not in SCs:
        SCs[perplexity]={}
    SCs[perplexity][learningRate]=[SC,nc]

# 2.2. build a matrix of SCs
perplexities=list(SCs.keys())
learningRates=list(SCs[perplexities[0]].keys())
perplexities.sort(); learningRates.sort()

Mlist=[]; Nlist=[]
for perplexity in perplexities:
    m=[]; n=[]
    for learningRate in learningRates:
        m.append(SCs[perplexity][learningRate][0]); n.append(SCs[perplexity][learningRate][1])
    Mlist.append(m); Nlist.append(n)
M=numpy.array(Mlist); N=numpy.array(Nlist)

# 2.3. build the figure
figureName='tsne.explorer.results.png'
figureTitle='Goodness of clustering after tSNE'

matplotlib.pyplot.imshow(M,interpolation='none',cmap='viridis')
cb=matplotlib.pyplot.colorbar(orientation='vertical',fraction=0.05) 
cb.set_label(label='Silhouette coefficient',size=16)
cb.ax.tick_params(labelsize=16)
matplotlib.pyplot.grid(False)

# setting the numbers
x=-0.1
y=0.1
deltax=1.
deltay=1.
for i in range(len(N)):
    for j in range(len(N[0])):
        value=str(N[i][j])
        matplotlib.pyplot.text(x+deltax*j,y+deltay*i,value,fontsize=10,color='black')

matplotlib.pyplot.xticks(range(len(learningRates)),learningRates,size=20,rotation=90)
matplotlib.pyplot.yticks(range(len(perplexities)),perplexities,size=20)

matplotlib.pyplot.xlabel('Learning rate',fontsize=18)
matplotlib.pyplot.ylabel('Perplexity',fontsize=18)
       
matplotlib.pyplot.tick_params(axis='x',which='both',bottom='off',top='off')
matplotlib.pyplot.tick_params(axis='y',which='both',right='off',left='off')
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.title(figureTitle,fontsize=18)
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(figureName)

matplotlib.pyplot.clf()

    
