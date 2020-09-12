if 'hard' in config['qc']['strategy']:
	label='hard'
else:
	label='outly'


if config['norm']['perSample']:
	if 'cca' in config['norm']['findAnchor']:
		label1='normPerSample_cca'
	else:
		label1='normPerSample_rpca'
else:
	label1='norm'


rule normalize:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+label+'flted_GeneUMIRatiogt0.8.rds')
	output:
		rds=path.join(config['dir']['data'],'{cancer}_{celltype}_'+label+'flted_'+label1+'.rds'),
		pca=path.join(config['dir']['plot'],'normalization','{cancer}_{celltype}_'+label+'flted_'+label1+'_pca.jpg'),
		tsne=path.join(config['dir']['plot'],'normalization','{cancer}_{celltype}_'+label+'flted_'+label1+'_tsne.jpg'),
		umap=path.join(config['dir']['plot'],'normalization','{cancer}_{celltype}_'+label+'flted_'+label1+'_umap.jpg'),
		combine=path.join(config['dir']['plot'],'normalization','{cancer}_{celltype}_'+label+'flted_'+label1+'_tsneVSumap.jpg')
	params:
		ccgenes=path.join(config['dir']['resources'],'cycle.rda'),
		sctPreNorm=config['norm']['sctPreNorm'],
		normPerSample=config['norm']['perSample'],
		nPC=config['norm']['maxPC'],
		reduction=config['norm']['findAnchor']		
	threads:
		10
	script:
		"../scripts/normalize.R"


rule normalize_all:
	input:
		expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+label+'flted_'+label1+'.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['plot'],'normalization','{cancer}_{celltype}_'+label+'flted_'+label1+'_pca.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['plot'],'normalization','{cancer}_{celltype}_'+label+'flted_'+label1+'_tsne.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['plot'],'normalization','{cancer}_{celltype}_'+label+'flted_'+label1+'_umap.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['plot'],'normalization','{cancer}_{celltype}_'+label+'flted_'+label1+'_tsneVSumap.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],'normalize_all.finish')
	shell:
		"touch {output}"



if config['norm']['sctPreNorm']: 
	label2='_sctPreNorm_'
else:
	label2='_'


rule checkCellCycle:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_'+label+'flted_GeneUMIRatiogt0.8.rds')
	output:
		split=path.join(config['dir']['plot'],'cellcycle','{cancer}_{celltype}_'+label+'flted'+label2+'cellcycle_split.jpg'),
		overlay=path.join(config['dir']['plot'],'cellcycle','{cancer}_{celltype}_'+label+'flted'+label2+'cellcycle_overlay.jpg')
	params:
		ccgenes=path.join(config['dir']['resources'],'cycle.rda'),
		sctPreNorm=config['norm']['sctPreNorm']
	script:
		"../scripts/check_cellcycle.R"



rule checkCellCycle_all:
	input:
		expand(path.join(config['dir']['plot'],'cellcycle','{cancer}_{celltype}_'+label+'flted'+label2+'cellcycle_split.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['plot'],'cellcycle','{cancer}_{celltype}_'+label+'flted'+label2+'cellcycle_overlay.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],'checkCellCycle_all.finish')		
	shell:
		"touch {output}"



