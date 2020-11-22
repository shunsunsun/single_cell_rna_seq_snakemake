#!/bin/bash

nCPU=$1
scriptName=$2  #scanorama_clust.R

source /gpfs2/gaog_pkuhpc/users/liny/tools/miniconda3/bin/activate gej_sc
export LD_LIBRARY_PATH=/home/gaog_pkuhpc/miniconda3/lib/:$LD_LIBRARY_PATH
source /appsnew/source/R-4.0.2share.sh

Rscript $scriptName $nCPU
