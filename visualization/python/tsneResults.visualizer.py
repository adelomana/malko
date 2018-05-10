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
for result in results:
    print(result)
