#if 'hard' in config['qc']['strategy']:
#	inData='hardflted'
#if 'outly' in config['qc']['strategy']:
#	inData='outlyflted'
#if 'pipe' in config['qc']['strategy']:
#	inData='pipeCompflted'


rule cca_integration:
	input:
		#path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_sctNorm.rds')
		path.join(config['dir']['data'],'{dataset}_QCed_sctNorm.rds')
	output:
		rds=path.join(config['dir']['data'],'{dataset}_QCed_sctNorm_BatchCCA.rds')
		#rds=path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_BatchS3.rds')
		#pca=path.join(config['dir']['plot'],'integrate','{cancer}_{celltype}_'+inData+'_BatchS3_pca.jpg'),
		#tsne=path.join(config['dir']['plot'],'integrate','{cancer}_{celltype}_'+inData+'_BatchS3_tsne.jpg'),
		#umap=path.join(config['dir']['plot'],'integrate','{cancer}_{celltype}_'+inData+'_BatchS3_umap.jpg'),
	params:
		#nPC=config['batchRm']['seurat3_nPC'],
		nFeature=config['batchRm']['nFeature'],
		anchor=config['batchRm']['seurat3_anchor']		
	threads:
		20
	script:
		"../scripts/step7_batchRm_seurat3.R"


#if 'sct' in config['batchRm']['norm4harmony']:
#	hmynorm='sctNorm'
#else:
#	hmynorm='logNormScale'


rule harmony:
	input:
		#path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_'+norm+'.rds')	
		#path.join(config['dir']['data'],'{dataset}_QCed_'+hmynorm+'.rds')
		path.join(config['dir']['data'],'{dataset}_QCed_sctNorm.rds')
	output:
		#rds=path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_'+norm+'_BatchHmy.rds')
		path.join(config['dir']['data'],'{dataset}_QCed_sctNorm_BatchHmy.rds')
	#params:
	#	#norm4harmony=hmynorm,
	#	theta=config['batchRm']['harmony_theta'],
	#	nclust=config['batchRm']['harmony_nclust'],
	#	max_it_clust=config['batchRm']['harmony_max_iter_cluster']
	threads:
		20
	script:
		"../scripts/step7_batchRm_harmony.R"


#if 'sct' in config['batchRm']['norm4fastMNN']:
#	fastMNNnorm='sctNorm'
#else:
#	fastMNNnorm='logNormScale'


rule fastMNN:
	input:
		#path.join(config['dir']['data'],'{dataset}_QCed_'+fastMNNnorm+'.rds')
		path.join(config['dir']['data'],'{dataset}_QCed_scNorm.rds')
	output:
		#path.join(config['dir']['data'],'{dataset}_QCed_'+fastMNNnorm+'_BatchfastMNN.rds')
                path.join(config['dir']['data'],'{dataset}_QCed_scNorm_BatchfastMNN.rds')
	params:
		#norm4fastMNN=fastMNNnorm,
		nfeatures=config['batchRm']['nFeature']
	threads:
		5
	script:
		"../scripts/step7_batchRm_fastMNN.R"


rule scanorama:
	input:
		path.join(config['dir']['data'],'{dataset}_QCed_sctNorm.rds')
	output:
		path.join(config['dir']['data'],'{dataset}_QCed_sctNorm_BatchScano.rds')
	threads:
		20
	script:
		"../scripts/step7_batchRm_scanorama.R"


rule combat:
	input:
		path.join(config['dir']['data'],'{dataset}_QCed_sctNorm.rds')
	output:
		path.join(config['dir']['data'],'{dataset}_QCed_sctNorm_BatchCombat.rds')
	#params:
	#	covariates="condition+cancerType"
	threads:
		20
	script:
		"../scripts/step7_batchRm_comBat.R"


def get_batch_input(wildcards):
	wanted = []
	if config['batchRm']['runSeurat3'] == True:
	#	wanted.extend(expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_BatchS3.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
		wanted.extend(expand(path.join(config['dir']['data'],'{dataset}_QCed_sctNorm_BatchS3.rds'),dataset=['ESCC','GEJ']))
	if config['batchRm']['runHarmony'] == True:
	#	wanted.extend(expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_'+hmynorm+'_BatchHmy.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
		wanted.extend(expand(path.join(config['dir']['data'],'{dataset}_QCed_'+hmynorm+'_BatchHmy.rds'),dataset=['ESCC','GEJ']))
	if config['batchRm']['runfastMNN'] == True:
		wanted.extend(expand(path.join(config['dir']['data'],'{dataset}_QCed_'+fastMNNnorm+'_BatchfastMNN.rds'),dataset=['ESCC','GEJ']))
	#if config['batchRm']['runCombat'] == True:
	#	wanted.extend(expand(path.join(config['dir']['data'],'{dataset}_QCed_sctNorm_BatchCombat.rds'),dataset=['ESCC','GEJ']))
	#if config['batchRm']['runScanorama'] == True:
	#	wanted.extend(expand(path.join(config['dir']['data'],'{dataset}_QCed_sctNorm_BatchScano.rds'),dataset=['ESCC','GEJ']))
	return wanted


rule batch_removal_all:
	input:
		get_batch_input
	output:
		path.join(config['dir']['log'],'batchRm_all.finish')
	shell:
		"touch {output}"


rule combat_all:
	input:
		expand(path.join(config['dir']['data'],'{dataset}_QCed_sctNorm_BatchCombat.rds'),dataset=['ESCC','GEJ'])
	output:
		path.join(config['dir']['log'],'combat.finish')
	shell:
		"touch {output}"


rule scanorama_all:
	input:
		expand(path.join(config['dir']['data'],'{dataset}_QCed_sctNorm_BatchScano.rds'),dataset=['ESCC','GEJ'])
	output:
		path.join(config['dir']['log'],'scanorama.finish')
	shell:
		"touch {output}"


rule cca_all:
	input:
		expand(path.join(config['dir']['data'],'{dataset}_QCed_sctNorm_BatchCCA.rds'),dataset=['ESCC','GEJ'])
	output:
		path.join(config['dir']['log'],'cca.finish')
	shell:
		"touch {output}"


rule harmony_all:
	input:
		expand(path.join(config['dir']['data'],'{dataset}_QCed_sctNorm_BatchHmy.rds'),dataset=['ESCC','GEJ'])
	output:
		path.join(config['dir']['log'],'harmony.finish')
	shell:
		"touch {output}"
