import os,time,numpy,datetime,pickle,sys
import multiprocessing,multiprocessing.pool
import matplotlib,matplotlib.pyplot

def launcher(iteration):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t iteration {}".format(iteration)))
    os.system('python worker.py -i {}'.format(iteration))

    return None

###
### MAIN
###

# 0. user defined variables
iterations=numpy.arange(80)
numberOfThreads=8

# 1. original analysis
#launcher(-1)

# 2. exploration
#print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t exploration"))
#hydra=multiprocessing.pool.Pool(numberOfThreads)
#hydra.map(launcher,iterations)

# 3. figure building
print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t figure building"))
figureFile='figure.pdf'

# read original
iteration=-1
jar='jars/jar.iter{}.pickle'.format(iteration)
f=open(jar,'rb')
results=pickle.load(f)
f.close()

neighbors=[]; originalK=[]; originalSS=[]; originalCHI=[]; originalVI=[]
for i in range(len(results)):
    if results[i][0] < 40:
        neighbors.append(results[i][0])
        originalK.append(results[i][1])
        originalSS.append(results[i][2])
        originalCHI.append(results[i][3])
        originalVI.append(results[i][4])
   
index=numpy.argmax(originalSS)
print('SS says nei = {}; K = {}'.format(neighbors[index],originalK[index]))
index=numpy.argmax(originalCHI)
print('CHI says nei = {}; K = {}'.format(neighbors[index],originalK[index]))
index=numpy.argmin(originalVI)
print('VI says nei = {}; K = {}'.format(neighbors[index],originalK[index]))

# read iterations
allK=[];allSS=[];allCHI=[];allVI=[]
for iteration in range(80):
    jar='jars/jar.iter{}.pickle'.format(iteration)
    f=open(jar,'rb')
    results=pickle.load(f)
    f.close()
    
    k=[]; SS=[]; CHI=[]; VI=[]
    for i in range(len(results)):
        if results[i][0] < 40:
            k.append(results[i][1])
            SS.append(results[i][2])
            CHI.append(results[i][3])
            VI.append(results[i][4])
    allK.append(k); allSS.append(SS); allCHI.append(CHI); allVI.append(VI)

meanK=numpy.mean(allK,axis=0)
meanSS=numpy.mean(allSS,axis=0)
meanCHI=numpy.mean(allCHI,axis=0)
meanVI=numpy.mean(allVI,axis=0)

fig,ax1=matplotlib.pyplot.subplots()
line1=ax1.plot(neighbors,originalK,'-',color='black')
#ax1.plot(neighbors,meanK,'-',color='black',lw=2)
#top=numpy.array(meanK)+numpy.array(stdK); bottom=numpy.array(meanK)-numpy.array(stdK)
#ax1.fill_between(neighbors,top,bottom,color='black',alpha=1/5,lw=0)

ax2=ax1.twinx()
ax2.set_yticks([])
line2=ax2.plot(neighbors,originalSS,'-',color='red',label='SS')
#ax2.plot(neighbors,meanSS,'-',color='red',lw=2)
#top=numpy.array(meanSS)+numpy.array(stdSS); bottom=numpy.array(meanSS)-numpy.array(stdSS)
#ax2.fill_between(neighbors,top,bottom,color='red',alpha=1/5,lw=0)

ax3=ax1.twinx()
ax3.set_yticks([])
line3=ax3.plot(neighbors,originalCHI,'-',color='blue',label='CHI')
#ax3.plot(neighbors,meanCHI,'-',color='blue',lw=2)
#top=numpy.array(meanCHI)+numpy.array(stdCHI); bottom=numpy.array(meanCHI)-numpy.array(stdCHI)
#ax3.fill_between(neighbors,top,bottom,color='blue',alpha=1/5,lw=0)

ax4=ax1.twinx()
ax4.set_yticks([])
line4=ax4.plot(neighbors,originalVI,'-',color='green',label='VI')
#ax4.plot(neighbors,meanVI,'-',color='green',lw=2)
#top=numpy.array(meanVI)+numpy.array(stdVI); bottom=numpy.array(meanVI)-numpy.array(stdVI)
#ax4.fill_between(neighbors,top,bottom,color='green',alpha=1/5,lw=0)

lns = line2+line3+line4
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs,loc=6)

#ax1.set_xlim([5,100])
ax1.set_ylabel('Partitions')
ax1.set_xlabel('Neighbors')
ax1.grid(linestyle=':',alpha=1/3)
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(figureFile)
matplotlib.pyplot.clf()
    
# 4. final message
#print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t completed"))
