#!/bin/bash

#$ -N tsneA
#$ -o /proj/omics4tb/alomana/scratch/messages.A.o.txt
#$ -e /proj/omics4tb/alomana/scratch/messages.A.e.txt
#$ -pe smp 40
#$ -S /bin/bash

cd /users/alomana
source .bash_profile

time python /proj/omics4tb/alomana/projects/mscni/src/deployment/tsne.explorer.deployment.A.py
