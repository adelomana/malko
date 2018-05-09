###
### This script visualizes the results from tsne explorer.
###

import sys,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

# 1. revover daata
f=open(similarityJar,'rb')
    similarities=pickle.load(f)
    f.close()

# 2. build figure
