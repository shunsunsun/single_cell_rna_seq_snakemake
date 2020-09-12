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

##### load rules #####

include: "rules/merge_cells.smk"
include: "rules/qc.smk"
include: "rules/filtration.smk"
#include: "rules/cell-cycle.smk"
include: "rules/normalization.smk"
#include: "rules/variance.smk"
#include: "rules/cell-type.smk"
#include: "rules/diffexp.smk"
