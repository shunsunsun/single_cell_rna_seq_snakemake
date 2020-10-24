import pandas as pd
from os import path

########## load config ############

configfile: "config/config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

tmp = pd.read_table(config['manifest']['gej_immun'],header=None)
GEJ_IM = list(tmp[0])
tmp = pd.read_table(config['manifest']['gej_epi'],header=None)
GEJ_EP = list(tmp[0])
tmp = pd.read_table(config['manifest']['escc_immun'],header=None)
ESCC_IM = list(tmp[0])
tmp = pd.read_table(config['manifest']['escc_epi'],header=None)
ESCC_EP = list(tmp[0])

ESCC_SAMPLES=list(set().union(ESCC_EP,ESCC_IM))
GEJ_SAMPLES=list(set().union(GEJ_EP,GEJ_IM))
SAMPLES=list(set().union(ESCC_SAMPLES,GEJ_SAMPLES))


##### load rules #####

include: "rules/merge_samples.smk"
include: "rules/qc.smk"
include: "rules/filtration.smk"
include: "rules/normalization.smk"
include: "rules/batchRm.smk"
include: "rules/dimensionReduct.smk"
include: "rules/clustering.smk"
#include: "rules/clust_qc.smk"
#include: "rules/variance.smk"
#include: "rules/cell-type.smk"
#include: "rules/diffexp.smk"
