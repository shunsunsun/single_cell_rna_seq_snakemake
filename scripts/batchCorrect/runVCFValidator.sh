#!/bin/bash

in=$1
out=$2

../tools/EBI_EVA_VCF_submission/vcf_validator_linux -i $in > $out
