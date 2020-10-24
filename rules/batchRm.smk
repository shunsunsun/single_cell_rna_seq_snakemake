if 'hard' in config['qc']['strategy']:
	inData='hardflted'
if 'outly' in config['qc']['strategy']:
	inData='outlyflted'
if 'pipe' in config['qc']['strategy']:
	inData='pipeCompflted'


rule seurat3_integration:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_sctNorm.rds')
	output:
		rds=path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_BatchS3.rds')
		#pca=path.join(config['dir']['plot'],'integrate','{cancer}_{celltype}_'+inData+'_BatchS3_pca.jpg'),
		#tsne=path.join(config['dir']['plot'],'integrate','{cancer}_{celltype}_'+inData+'_BatchS3_tsne.jpg'),
		#umap=path.join(config['dir']['plot'],'integrate','{cancer}_{celltype}_'+inData+'_BatchS3_umap.jpg'),
	params:
		#nPC=config['batchRm']['seurat3_nPC'],
		nFeature=config['batchRm']['seurat3_nFeature'],
		anchor=config['batchRm']['seurat3_anchor']		
	threads:
		20
	script:
		"../scripts/step7_batchRm_seurat3.R"


if 'sct' in config['batchRm']['norm4harmony']:
	hmynorm='sctNorm'
else:
	hmynorm='logNormScale'


rule harmony:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_'+hmynorm+'.rds')	
	output:
		rds=path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_'+hmynorm+'_BatchHmy.rds')
	params:
		norm4harmony=hmynorm,
		theta=config['batchRm']['harmony_theta'],
		nclust=config['batchRm']['harmony_nclust'],
		max_it_clust=config['batchRm']['harmony_max_iter_cluster']
	threads:
		20
	script:
		"../scripts/step7_batchRm_harmony.R"


def get_batch_input(wildcards):
	wanted = []
	if config['batchRm']['runSeurat3'] == True:
		wanted.extend(expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_BatchS3.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['batchRm']['runHarmony'] == True:
		wanted.extend(expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_'+hmynorm+'_BatchHmy.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	return wanted


rule batch_removal_all:
	input:
		get_batch_input
	output:
		path.join(config['dir']['log'],'batchRm_all.finish')
	shell:
		"touch {output}"


