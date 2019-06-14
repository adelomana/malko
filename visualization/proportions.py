import numpy,sys
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

# read data
proportions={}

f=open('worldColors.eigengenes.txt','r')
for line in f:
    v=line.split('\t')
    dayLabel=v[0]
    elements=v[1].split(',')
    states=[int(element) for element in elements]
    proportions[dayLabel]=states
f.close()

print(proportions)


# compute proportions
cellStates=list(range(0,11))


dynamics=[]
for dayLabel in proportions.keys():
    print(dayLabel)
    states=proportions[dayLabel]
    print(states[:50])
    
    numberOfCells=len(states)
    print(numberOfCells)


    percentages=[]
    for element in cellStates:
        p=int((100*states.count(element))/numberOfCells)
        percentages.append(p)
    
    dynamics.append(percentages)
    
    print()

# plot
D=numpy.array(dynamics)
D=numpy.transpose(D)

import scipy,scipy.cluster,scipy.cluster.hierarchy
scipy.cluster.hierarchy.linkage(D)

sys.exit()

matplotlib.pyplot.imshow(D,cmap='viridis')
cb=matplotlib.pyplot.colorbar(label='Cell state %',orientation='vertical',fraction=0.1) # need to improve label font size
cb.ax.tick_params(labelsize=20)
cb.ax.set_yticklabels([2**element for element in [3,4,5,6]])
matplotlib.pyplot.grid(False)

    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.tick_params(axis='x',which='both',bottom='off',top='off')
matplotlib.pyplot.tick_params(axis='y',which='both',right='off',left='off')
matplotlib.pyplot.axes().set_aspect('equal')
#matplotlib.pyplot.xlabels(['control',d])

labels=['ctl','d3','d6','d13','d17','d24']

matplotlib.pyplot.xticks(range(len(labels)),labels,size=18,rotation=33)
matplotlib.pyplot.ylabel('Cell states')

matplotlib.pyplot.savefig('dynamics.pdf')
print(dynamics)


