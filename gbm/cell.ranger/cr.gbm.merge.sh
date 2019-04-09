#!/bin/bash

#$ -N mergeGBM
#$ -o /proj/omics4tb2/alomana/scratch/messages.GBM.merge.o.txt
#$ -e /proj/omics4tb2/alomana/scratch/messages.GBM.merge.e.txt
#$ -pe smp 40
#$ -l hostname=baliga2
#$ -S /bin/bash

cd /users/alomana
source .bash_profile 

echo "Original paths given my bash_profile..."
echo $PATH
echo ""

echo "Adding extra paths for Cell Ranger..."
source /proj/omics4tb2/alomana/software/cellRanger/cellranger-3.0.2/sourceme.bash
echo $PATH
echo ""

echo "Adding path to bcl2fastq..."
export PATH="/usr/local/bcl2fastq-v2.20.0.422/bin:$PATH"
echo $PATH
echo ""

echo "where am i?"
uname -a
echo ""

echo "Running count..." 
cd /proj/omics4tb2/alomana/projects/gbm/results/
time cellranger aggr --id=aggregated --csv=/proj/omics4tb2/alomana/projects/gbm/data/aggregation.csv --normalize=mapped --localcores=40 --localmem=90
echo "... run completed."
