!/bin/bash

gliny2
cd temp
for f in /gpfs2/gaog_pkuhpc/users/liny/strelkaResult/T[1-9]*/strelka.whole.flt.vcf.gz; do sample=$(basename $(dirname $f)); cp $f ${sample}.vcf.gz; done &
for f in T*vcf.gz; do gunzip $f; done &
for f in *vcf; do sample=${f/.vcf/}; normal=${sample/T/N}; sed -i "s/TUMOR/${sample}/g" $f; sed -i "s/NORMAL/${normal}/g" $f; done &
source /gpfs2/gaog_pkuhpc/users/liny/tools/miniconda3/bin/activate hcc_gg
for f in ./T*vcf; do bgzip $f; tabix ${f}.gz; done &
for f in ./T*vcf.gz; do sbatch -p cn-long -A gaog_g1 --qos=gaogcnl -N 1 -n 1 -c 1 ./runVCFValidator.sh $f ${f}.out; sleep 1; done
for o in ./T*summary*txt; do msg=$(cat $o | grep "the input file is valid"); if [[ -z $msg ]]; then echo $o; fi; done
for f in T*vcf.gz; do md5sum $f | awk -F' ' -vOFS='\t' '{print $2,"vcf",$1}' >> result; md5sum ${f}.tbi | awk -F' ' -vOFS='\t' '{print $2,"tabix",$1}' >> result; done &

