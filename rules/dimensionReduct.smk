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


rule pre_biscuit:
	input:
		path.join(config['dir']['data'],'{dataset}_QCed.rds')
	output:
		path.join(config['dir']['scripts'],'BISCUIT','{dataset}_QCed_4biscuit.rda')
	script:
		"../scripts/BISCUIT/prepare_BISCUIT_input.R"


rule biscuit:
	input:
		path.join(config['dir']['scripts'],'BISCUIT','{dataset}_QCed_4biscuit.rda')
	output:
		path.join(config['dir']['scripts'],'BISCUIT','{dataset}','plots/extras/Output_values.txt')
	params:
		infile='{dataset}_QCed_4biscuit.rda',
		outdir='{dataset}',
		alpha=0.5
	threads:
		30
	script:
		"../scripts/BISCUIT/start_file.R"


rule biscuit_all:
	input:
		expand(path.join(config['dir']['scripts'],'BISCUIT','{dataset}','plots/extras/Output_values.txt'),dataset=['ESCC','GEJ'])
	output:
		path.join(config['dir']['log'],'biscuit.finish')
	script:
		"touch {output}"


rule pca:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_'+norm+'.rds')
	output:
		rds=path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_'+norm+'_pca.rds'),
		optimPC=path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_'+norm+'_optimalPC.txt')
	params:
		plot=path.join(config['dir']['plot'],'dimenReduct'),
		nhvg=config['batchRm']['nFeature']				
	script:
		"../scripts/step8_dimRed_vstpca.R"


#ndim=str(config['glmpca']['ndim'])


rule glm_pca:
	input:
		#path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'.rds')		
		#path.join(config['dir']['data'],'{dataset}_QCed.rds')
		#path.join(config['dir']['data'],'GEJ_QCed_sctNorm_BatchCCA_clustStab','{celltype}.rds')
		path.join(config['dir']['data'],'ESCC_QCed_sctNorm_BatchHmy_clustStab','{celltype}.rds')
	output:
		#path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_glmpca_ndim'+ndim+'.rds')
		#path.join(config['dir']['data'],'{dataset}_QCed_glmpcaReduct.rds')
		#path.join(config['dir']['data'],'GEJ_QCed_sctNorm_BatchCCA_clustStab','{celltype}_glmpca.rds')
		path.join(config['dir']['data'],'ESCC_QCed_sctNorm_BatchHmy_clustStab','{celltype}_glmpca.rds')
	params:
		plot=path.join(config['dir']['plot'],'dimenReduct'),
		nhvg=config['batchRm']['nFeature'],
		ndims="10_20"
		#ndim=config['glmpca']['ndim']
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
		#expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+flt+'_glmpca_ndim'+ndim+'.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
		#expand(path.join(config['dir']['data'],'{dataset}_QCed_glmpcaReduct.rds'),dataset=['ESCC','GEJ'])
                #expand(path.join(config['dir']['data'],'GEJ_QCed_sctNorm_BatchCCA_clustStab','{celltype}_glmpca.rds'),celltype=['Monocytes','TCells','BCells','ImmuneCells','StromalCells','EpithelialCells'])
                expand(path.join(config['dir']['data'],'ESCC_QCed_sctNorm_BatchHmy_clustStab','{celltype}_glmpca.rds'),celltype=['Monocytes','TCells','BCells','ImmuneCells','StromalCells'])
	output:
		path.join(config['dir']['log'],'glmpca_all.finish')
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
	

