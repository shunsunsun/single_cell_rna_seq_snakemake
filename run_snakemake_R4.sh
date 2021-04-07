#!/bin/bash

task=$1
keepGoing=$2
cpu=$3
njob=$4
cluster=$5

##Run the workflow at the directory that hosts the Snakefile
##sbatch -p cn-long -A gaog_g1 --qos=gaogcnl -N 1 -n 1 -c 1 ./run_snakemake.sh \
##<rule_name> <whether-keep-going> <assigned cpus> <number of jobs> <partition>
##Example: sbatch -p cn-icg -A gaog_g1 --qos=gaogcnicg -N 1 -n 1 -c 1 snakemake.sh all n 5 50 icg

source /gpfs2/gaog_pkuhpc/users/liny/tools/miniconda3/bin/activate gej_sc
export LD_LIBRARY_PATH=/gpfs2/gaog_pkuhpc/users/liny/tools/miniconda3/envs/gej_sc/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/gaog_pkuhpc/miniconda3/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/appsnew/usr/python/Miniconda/miniconda3-4.7.10/envs/gdal/lib/:$LD_LIBRARY_PATH
source /appsnew/source/R-4.0.2share.sh
#cellphone DB
#source /gpfs2/gaog_pkuhpc/users/liny/tools/miniconda3/bin/activate cellphonedb

if [[ ${keepGoing} == "y" ]]; then
	option1="--keep-going"
fi
#if [[ ${cluster} == "long" ]]; then
	#sbcmd="sbatch -p cn-long -A gaog_g1 --qos=gaogcnl -N 1 -n 1 -c ${cpu}"
if [[ ${cluster} == "bio" ]]; then
	sbcmd="sbatch -p cn_bio -A gaog_cg6 --qos=gaogbioq -N 1 -n 1 -c ${cpu}"
elif [[ ${cluster} == "icg" ]]; then
	sbcmd="sbatch -p cn_icg -A gaog_g1 --qos=gaogcnicg -N 1 -n 1 -c ${cpu}"
elif [[ ${cluster} == "fati" ]]; then
        sbcmd="sbatch -p fat_icg -A gaog_g1 --qos=gaogfaticg -N 1 -n 1 -c ${cpu}"
fi
snakemake ${task} ${option1} --max-jobs-per-second 1 --jobs ${njob} --cluster "${sbcmd}" --latency-wait 120
