import pickle,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

figureFile='figures/figure.003.pdf'
jar='jars/exploration.results.003.pickle'
f=open(jar,'rb')
results=pickle.load(f)
f.close()

neighbors=[]
originalK=[]; originalSS=[]; originalCHI=[]; originalVI=[]
meanK=[]; meanSS=[]; meanCHI=[]; meanVI=[]
stdK=[]; stdSS=[]; stdCHI=[]; stdVI=[]

for iteration in results:

    original=iteration[0]
    meanResult=iteration[1]
    stdResult=iteration[2]

    neighbors.append(original[0])

    originalK.append(original[1]); originalSS.append(original[2]); originalCHI.append(original[3]); originalVI.append(original[4])
    meanK.append(meanResult[1]); meanSS.append(meanResult[2]); meanCHI.append(meanResult[3]); meanVI.append(meanResult[4])
    stdK.append(stdResult[1]); stdSS.append(stdResult[2]); stdCHI.append(stdResult[3]); stdVI.append(stdResult[4])    

fig,ax1=matplotlib.pyplot.subplots()
line1=ax1.plot(neighbors,originalK,'o',color='black',mew=0)
ax1.plot(neighbors,meanK,'-',color='black',lw=2)
top=numpy.array(meanK)+numpy.array(stdK); bottom=numpy.array(meanK)-numpy.array(stdK)
ax1.fill_between(neighbors,top,bottom,color='black',alpha=1/5,lw=0)

ax2=ax1.twinx()
ax2.set_yticks([])
line2=ax2.plot(neighbors,originalSS,'o',color='red',label='SS',mew=0)
ax2.plot(neighbors,meanSS,'-',color='red',lw=2)
top=numpy.array(meanSS)+numpy.array(stdSS); bottom=numpy.array(meanSS)-numpy.array(stdSS)
ax2.fill_between(neighbors,top,bottom,color='red',alpha=1/5,lw=0)

ax3=ax1.twinx()
ax3.set_yticks([])
line3=ax3.plot(neighbors,originalCHI,'o',color='blue',label='CHI',mew=0)
ax3.plot(neighbors,meanCHI,'-',color='blue',lw=2)
top=numpy.array(meanCHI)+numpy.array(stdCHI); bottom=numpy.array(meanCHI)-numpy.array(stdCHI)
ax3.fill_between(neighbors,top,bottom,color='blue',alpha=1/5,lw=0)

ax4=ax1.twinx()
ax4.set_yticks([])
line4=ax4.plot(neighbors,originalVI,'o',color='green',label='VI',mew=0)
ax4.plot(neighbors,meanVI,'-',color='green',lw=2)
top=numpy.array(meanVI)+numpy.array(stdVI); bottom=numpy.array(meanVI)-numpy.array(stdVI)
ax4.fill_between(neighbors,top,bottom,color='green',alpha=1/5,lw=0)

ax1.set_ylabel('Partitions')
ax1.set_xlabel('Neighbors')
ax1.grid(linestyle=':',alpha=1/3)
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(figureFile)
matplotlib.pyplot.clf()


