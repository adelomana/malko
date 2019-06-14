#!/bin/bash

#cp /Volumes/omics4tb2/alomana/projects/mscni/data/original/output\ of\ cell\ ranger/*.zip .

unzip 20765_M397_control.zip
#unzip 21729_day3.zip
#unzip 21766_day6.zip
#unzip 22077_day13.zip
#unzip 22152_day17.zip
#unzip 18324_day24.zip

time velocyto run10x -m /Users/alomana/scratch/hsa_hg19_rmsk.gtf /Users/alomana/scratch/20765_M397_control /Users/alomana/scratch/genes.gtf --samtools-memory 1500 --samtools-threads 8 -v > messages.control.txt

time velocyto run10x -m /Users/alomana/scratch/hsa_hg19_rmsk.gtf /Users/alomana/scratch/21729_day3 /Users/alomana/scratch/genes.gtf --samtools-memory 1500 --samtools-threads 8 -v > messages.day.3.txt

time velocyto run10x -m /Users/alomana/scratch/hsa_hg19_rmsk.gtf /Users/alomana/scratch/21766_day6 /Users/alomana/scratch/genes.gtf --samtools-memory 1500 --samtools-threads 8 -v > messages.day.6.txt

time velocyto run10x -m /Users/alomana/scratch/hsa_hg19_rmsk.gtf /Users/alomana/scratch/22077_day13 /Users/alomana/scratch/genes.gtf --samtools-memory 1500 --samtools-threads 8 -v > messages.day.13.txt

time velocyto run10x -m /Users/alomana/scratch/hsa_hg19_rmsk.gtf /Users/alomana/scratch/22152_day17 /Users/alomana/scratch/genes.gtf --samtools-memory 1500 --samtools-threads 8 -v > messages.day.17.txt

time velocyto run10x -m /Users/alomana/scratch/hsa_hg19_rmsk.gtf /Users/alomana/scratch/18324_Yapeng_single_cell /Users/alomana/scratch/genes.gtf --samtools-memory 1500 --samtools-threads 8 -v > messages.day.24.txt
