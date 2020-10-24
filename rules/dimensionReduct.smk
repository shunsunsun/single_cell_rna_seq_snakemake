if 'hard' in config['qc']['strategy']:
        flt='hardflted'
if 'outly' in config['qc']['strategy']:
        flt='outlyflted'
if 'pipe' in config['qc']['strategy']:
        flt='pipeCompflted'
if 'sct' in config['norm']['method']:
        norm='sctNorm'
if 'log' in config['norm']['method']:
        norm='logNormScale'

rule pca:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_'+norm+'.rds')
	output:
		rds=path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_'+norm+'_pca.rds'),
		optimPC=path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_'+norm+'_optimalPC.txt')
	params:
		plot=path.join(config['dir']['plot'],'dimenReduct'),
		nhvg=config['featSelect']['nFeature']				
	script:
		"../scripts/step8_dimRed_vstpca.R"


ndim=str(config['glmpca']['ndim'])


rule glm_pca:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'.rds')		
	output:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_glmpca_ndim'+ndim+'.rds'),
	params:
		plot=path.join(config['dir']['plot'],'dimenReduct'),
		nhvg=config['featSelect']['nFeature'],
		ndim=ndim
	threads:
		10
	script:
		"../scripts/step8_dimRed_devGLMpca.R"


rule pca_all:
	input:
		expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_'+norm+'_pca.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_'+norm+'_optimalPC.txt'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],'pca_all.finish')
	shell:
		"touch {output}"
		

rule glm_pca_all:
	input:
		expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_glmpca_ndim'+ndim+'.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],'glm_pca_all.finish')
	shell:
		"touch {output}"


rule select_pca:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_sctNorm_pca.rds')
	output:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_dimReduc.rds')
	params:
		glmpca=path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_glmpca'),
		pca_dr="pca_2000hvg",
		glm_dr="glmpca_dev2000"		
	script:
		"../scripts/step8b_selectDimReduc.R"


rule select_pca_all:
	input:
		expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_dimReduc.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],'dimRec_select.finish')
	shell:
		"touch {output}"
	

