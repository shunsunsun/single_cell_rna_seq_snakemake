if 'outly' in config['qc']['strategy']:
        flt='outlyflted'
if 'pipe' in config['qc']['strategy']:
        flt='pipeCompflted'


rule clust_seurat3:
        input:
                path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_dimReduc.rds')
        output:
                path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_seClust.rds')
        params:
                config['cluster']['nPC']
	threads:
		20
	script:
                "../scripts/step9_clust_seurat3.R"


rule clust_scran:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_dimReduc.rds')
	output:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_walkClust.rda')
	params:
		config['cluster']['nPC']
	threads:
		20
	script:
		"../scripts/step9_clust_scran.R"


def get_clust_res(wildcards):
	wanted = []
	if 'pipe' in config['qc']['strategy']:
		wanted.extend(expand(path.join(config['dir']['data'],'{cancer}_{celltype}_pipeCompflted_{method}.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM'],method=['seClust']))
		wanted.extend(expand(path.join(config['dir']['data'],'{cancer}_{celltype}_pipeCompflted_{method}.rda'),cancer=['ESCC','GEJ'],celltype=['EP','IM'],method=['walkClust']))
	if 'outly' in config['qc']['strategy']:
		wanted.extend(expand(path.join(config['dir']['data'],'{group}_outlyflted_{method}.rds'),group=['ESCC_EP','ESCC_IM','GEJ_EP'],method=['seClust']))
		wanted.extend(expand(path.join(config['dir']['data'],'{group}_outlyflted_{method}.rda'),group=['ESCC_EP','ESCC_IM','GEJ_EP'],method=['walkClust']))
	return wanted


rule clust_all:
	input:
		get_clust_res
	output:
		path.join(config['dir']['log'],'clust.finish')
	shell:
		"touch {output}"


rule clust_diagnosis_seurat3:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_seClust.rds')
	output:
		path.join(config['dir']['log'],'{cancer}_{celltype}_'+flt+'_seClustEval.done')
	params:
		path.join(config['dir']['plot'],'clustering')
	script:
		"../scripts/step9_clustDiagnosis.R"


def get_clust_eval(wildcards):
	wanted = []
	if 'pipe' in config['qc']['strategy']:
		wanted.extend(expand(path.join(config['dir']['log'],'{cancer}_{celltype}_pipeCompflted_{method}Eval.done'),cancer=['ESCC','GEJ'],celltype=['EP','IM'],method=['seClust','walkClust']))
		#wanted.extend(expand(path.join(config['dir']['log'],'{cancer}_{celltype}_pipeCompflted_{method}Eval.done'),cancer=['ESCC','GEJ'],celltype=['EP','IM'],method=['seClust']))
	if 'outly' in config['qc']['strategy']:
		wanted.extend(expand(path.join(config['dir']['log'],'{group}_outlyflted_{method}Eval.done'),group=['ESCC_EP','ESCC_IM','GEJ_EP'],method=['seClust','walkClust']))
		#wanted.extend(expand(path.join(config['dir']['log'],'{group}_outlyflted_{method}Eval.done'),group=['ESCC_EP','ESCC_IM','GEJ_EP'],method=['seClust']))
	return wanted


rule clust_diagnosis_scran:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_walkClust.rda')
	output:
		path.join(config['dir']['log'],'{cancer}_{celltype}_'+flt+'_walkClustEval.done')
	params:
		path.join(config['dir']['plot'],'clustering')
	script:
		"../scripts/step9_clustDiagnosis.R"


rule clust_diagnosis_all:
	input:
		get_clust_eval
	output:
		path.join(config['dir']['log'],flt+'_clust_eval.finish')
	shell:
		"touch {output}"


select_clust={'GEJ_EP':'walktrap_glmpca10_50','GEJ_IM':'walktrap_glmpca10_50','ESCC_IM':'seurat_glmpca20_0.6','ESCC_EP':'walktrap_glmpca20_30'}

rule select_clust_diagn:
	input:
		path.join(config['dir']['data'], '{tissue_cell_type}_pipeCompflted_seClust.rds')
	output:
		path.join(config['dir']['log'],'{tissue_cell_type}_pipeCompflted_stableClust.rpt')
	params:
		outPlot=path.join(config['dir']['plot'],'clustering'),
		clust=lambda wildcards: select_clust[wildcards.tissue_cell_type]
	script:
		"../scripts/step9_selectClustDiagnosis.R"


rule select_clust_diagn_all:
	input:
		expand(path.join(config['dir']['log'],'{tissue_cell_type}_pipeCompflted_stableClust.rpt'),tissue_cell_type=['ESCC_EP','GEJ_EP','ESCC_IM','GEJ_IM'])
	output:
		path.join(config['dir']['log'],'select_clust_diagn.finish')
	shell:
		"touch {output}"


rule clean_clust_res:
	input:
		path.join(config['dir']['data'], '{tissue_cell_type}_pipeCompflted_seClust.rds')
	output:
		path.join(config['dir']['data'], '{tissue_cell_type}_pipeCompflted_clustered.rds')
	params:
		clust=lambda wildcards: select_clust[wildcards.tissue_cell_type]
	script:
		"../scripts/step9b_cleanAfterClust.R"


rule clean_clust_all:
	input:
		expand(path.join(config['dir']['data'], '{tissue_cell_type}_pipeCompflted_clustered.rds'),tissue_cell_type=['ESCC_EP','ESCC_IM','GEJ_EP','GEJ_IM'])
		#expand(path.join(config['dir']['data'], '{tissue_cell_type}_pipeCompflted_clustered.rds'),tissue_cell_type=['ESCC_IM','GEJ_EP','GEJ_IM'])
	output:
		path.join(config['dir']['log'],'clean_clust.finish')
	shell:
		"touch {output}"
