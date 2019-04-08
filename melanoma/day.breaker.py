###
### This script separates original data into different files per day
###

import sys

# 0. user defined variables
originalDataFile='/Volumes/omics4tb2/alomana/projects/mscni/results/seurat.cleaned.data.csv'
outputDirectory='/Volumes/omics4tb2/alomana/projects/mscni/results/broken/'

# 1. read data
expression={}
days=[]
with open(originalDataFile,'r') as f:
    header=f.readline()
    v=header.split(',')
    elements=v[1:]
    repeatedDays=[element.split('_')[1] for element in elements]
    repeatedDays=[element.replace('\n','') for element in repeatedDays]

    uniqueDayLabels=list(set(repeatedDays))
    print(uniqueDayLabels)
    sys.exit()
    

    for line in f:
        v=line.split(',')
        geneName=v[0]
        expressionValues=v[1:]
        values=[float(element) for element in expressionValues]
        
        #print(v,len(v),v[:4])
        print(geneName,expressionValues[:10])
        for i in range(len(values)):
            dayLabel=repeatedDays[i]
            value=values[i]

            if dayLabel not in expression:
                expression[dayLabel]=[]

            

            print(len(repeatedDays),len(values))
            
    sys.exit()


# break it into smaller variables

# run tSNE

# save smaller files

# 2. write new files, splitted data per day
