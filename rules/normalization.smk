if 'hard' in config['qc']['strategy']:
        flt='hardflted'
if 'outly' in config['qc']['strategy']:
        flt='outlyflted'
if 'pipe' in config['qc']['strategy']:
	flt='pipeCompflted'


rule sctNormalize:
	input:
		#path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'.rds')
		path.join(config['dir']['data'],'GEJ_QCed_sctNorm_BatchCCA_clustStab','{celltype}.rds')
		#path.join(config['dir']['data'],'ESCC_QCed_sctNorm_BatchHmy_clustStab','{celltype}.rds')
	output:
		#path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_sctNorm.rds')
		path.join(config['dir']['data'],'GEJ_QCed_sctNorm_BatchCCA_clustStab','{celltype}_sctNorm.rds')
		#path.join(config['dir']['data'],'ESCC_QCed_sctNorm_BatchHmy_clustStab','{celltype}_sctNorm.rds')
	params:
		regressNum=config['norm']['regressNum'],		
		ccgenes=path.join(config['dir']['resources'],'cycle.rda'),
		sctPreNorm=config['norm']['sctPreNorm']
	threads:
		20
	script:
		#"../scripts/step6_scTransform.R"
		"../scripts/step6_scTransform_subset.R"


rule logNormalize:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'.rds')
	output:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_logNormScale.rds')
	params:
		regressNum=config['norm']['regressNum'],
		ccgenes=path.join(config['dir']['resources'],'cycle.rda')
	threads:
		20
	script:
		"../scripts/step6_logNormScale.R"


if 'sct' in config['norm']['method']:
	norm='sctNorm'
if 'log' in config['norm']['method']:
	norm='logNormScale'


rule whole_normalize_all:
	input:
		#expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_'+norm+'.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
		expand(path.join(config['dir']['data'],'GEJ_QCed_sctNorm_BatchCCA_clustStab','{celltype}_sctNorm.rds'),celltype=['Monocytes','TCells','BCells','ImmuneCells','StromalCells','EpithelialCells'])
		#expand(path.join(config['dir']['data'],'ESCC_QCed_sctNorm_BatchHmy_clustStab','{celltype}_sctNorm.rds'),celltype=['Monocytes','TCells','BCells','ImmuneCells','StromalCells','EpithelialCells'])
	output:
		path.join(config['dir']['log'],'norm_all.finish')
	shell:
		"touch {output}"


rule sctNorm_mergeData:
	input:
		path.join(config['dir']['data'],'{dataset}_QCed.rds')
	output:
		path.join(config['dir']['data'],'{dataset}_QCed_sctNorm.rds')
	params:
		regressNum=config['norm']['regressNum'],
		ccgenes=path.join(config['dir']['resources'],'cycle.rda'),
		sctPreNorm=config['norm']['sctPreNorm']
	threads:
		40
	script:
		"../scripts/step6_scTransform.R"


rule normalize_mergeData:
	input:
		expand(path.join(config['dir']['data'],'{dataset}_QCed_sctNorm.rds'),dataset=['ESCC','GEJ','IM'])
	output:
		path.join(config['dir']['log'],'norm_merged.finish')
	shell:
		"touch {output}"


rule split_sctNorm_ESCC_EP:
	input:
		path.join(config['dir']['data'],'ESCC_EP_'+flt+'.rds')
	output:
		path.join(config['dir']['log'],'ESCC_EP_'+flt+'_sctNorm.rpt'),
		expand(path.join(config['dir']['data'],'sctNormSeu_'+flt,'{sample}.sctNormSeu_'+flt+'.rds'),sample=ESCC_EP)
	params:
		regressNum=config['norm']['regressNum'],
		dir=path.join(config['dir']['data'],'sctNormSeuObj_'+flt),
		ccgenes=path.join(config['dir']['resources'],'cycle.rda'),
		sctPreNorm=config['norm']['sctPreNorm']
	threads:
		20
	script:
		"../scripts/step6_split_norm.R"


rule split_sctNorm_ESCC_IM:
	input:
		path.join(config['dir']['data'],'ESCC_IM_'+flt+'.rds')
	output:
		path.join(config['dir']['log'],'ESCC_IM_'+flt+'_sctNorm.rpt'),
		expand(path.join(config['dir']['data'],'sctNormSeu_'+flt,'{sample}.sctNormSeu_'+flt+'.rds'),sample=ESCC_IM)
	params:
		regressNum=config['norm']['regressNum'],
		dir=path.join(config['dir']['data'],'sctNormSeuObj_'+flt),
		ccgenes=path.join(config['dir']['resources'],'cycle.rda'),
		sctPreNorm=config['norm']['sctPreNorm']
	threads:
		20
	script:
		"../scripts/step6_split_norm.R"


rule split_sctNorm_GEJ_IM:
	input:
		path.join(config['dir']['data'],'GEJ_IM_'+flt+'.rds')
	output:
		path.join(config['dir']['log'],'GEJ_IM_'+flt+'_sctNorm.rpt'),
		expand(path.join(config['dir']['data'],'sctNormSeu_'+flt,'{sample}.sctNormSeu_'+flt+'.rds'),sample=GEJ_IM)
	params:
		regressNum=config['norm']['regressNum'],
		dir=path.join(config['dir']['data'],'sctNormSeuObj_'+flt),
		ccgenes=path.join(config['dir']['resources'],'cycle.rda'),
		sctPreNorm=config['norm']['sctPreNorm']
	threads:
		20
	script:
		"../scripts/step6_split_norm.R"


rule split_sctNorm_GEJ_EP:
	input:
		path.join(config['dir']['data'],'GEJ_EP_'+flt+'.rds')
	output:
		path.join(config['dir']['log'],'GEJ_EP_'+flt+'_sctNorm.rpt'),
		expand(path.join(config['dir']['data'],'sctNormSeu_'+flt,'{sample}.sctNormSeu_'+flt+'.rds'),sample=GEJ_EP)
	params:
		regressNum=config['norm']['regressNum'],
		dir=path.join(config['dir']['data'],'sctNormSeuObj_'+flt),
		ccgenes=path.join(config['dir']['resources'],'cycle.rda'),
		sctPreNorm=config['norm']['sctPreNorm']
	threads:
		20
	script:
		"../scripts/step6_split_norm.R"



##logNorm is conducted per cell, whereas sctNorm the entire dataset, so for the former, we can just normalize the whole instead of split --> norm --> merge
rule split_sctNorm_all:
	input:
		expand(path.join(config['dir']['log'],'{cancer}_{celltype}_'+flt+'_sctNorm.rpt'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['data'],'sctNormSeu_'+flt,'{sample}.sctNormSeu_'+flt+'.rds'),sample=SAMPLES)		
	output:
		path.join(config['dir']['log'],'split_sctNorm_all.finish')
	shell:
		"touch {output}"
