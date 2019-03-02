###
### This script separates original data into different files per day
###

import sys

# 0. user defined variables
originalDataFile='/Volumes/omics4tb2/alomana/projects/mscni/data/for.seurat/count.file.alldays.csv'
outputDirectory='/Volumes/omics4tb2/alomana/projects/mscni/data/per.day/'

# 1. read data
expression={}
days=[]
with open(originalDataFile,'r') as f:
    header=f.readline()
    v=header.split(',')
    elements=v[1:]
    repeatedDays=[element.split('_')[1] for element in elements]
    print(repeatedDays[:10])
    sys.exit()

    for line in f:
        v=line.split(',')
        print(v,len(v),v[:4])
        

# 2. write new files, splitted data per day
