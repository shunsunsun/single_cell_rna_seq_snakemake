if 'hard' in config['qc']['strategy']:
	label='hard'
else:
	label='outly'


rule plot_sctNorm:
	input:
		path.join(config['dir']['data'],'sctNormSeuObj','{sample}.sctNormSeu.rds')
	output:
		#pca=path.join(config['dir']['plot'],'sctNorm','{sample}.sctNorm.clust.pca.jpg'),
		#tsne=path.join(config['dir']['plot'],'sctNorm','{sample}.sctNorm.clust.tsne.jpg'),
		#umap=path.join(config['dir']['plot'],'sctNorm','{sample}.sctNorm.clust.umap.jpg'),
		qc=path.join(config['dir']['plot'],'sctNorm','{sample}.sctNorm.clust.qc.jpg')
	params:
		nPC=config['norm']['maxPC']
	threads:
		4
	script:
		"../scripts/clust_qc.R"


rule plot_rmLCsctNorm:
	input:
		path.join(config['dir']['data'],'rmLCsctNormSeuObj','{sample}.rmLCsctNormSeu.rds')
	output:
		#pca=path.join(config['dir']['plot'],'rmLCsctNorm','{sample}.rmLCsctNorm.clust.pca.jpg'),
		#tsne=path.join(config['dir']['plot'],'rmLCsctNorm','{sample}.rmLCsctNorm.clust.tsne.jpg'),
		#umap=path.join(config['dir']['plot'],'rmLCsctNorm','{sample}.rmLCsctNorm.clust.umap.jpg'),
		qc=path.join(config['dir']['plot'],'rmLCsctNorm','{sample}.rmLCsctNorm.clust.qc.jpg')		
	params:
		nPC=config['norm']['maxPC']
	threads:
		4
	script:
		"../scripts/clust_qc.R"


rule norm_plot_GEJ_EP:
	input:
		#expand(path.join(config['dir']['plot'],'sctNorm','{sample}.sctNorm.clust.{type}.jpg'),sample=GEJ_EP,type=['pca','tsne','umap','qc']),
		#expand(path.join(config['dir']['plot'],'rmLCsctNorm','{sample}.rmLCsctNorm.clust.{type}.jpg'),sample=GEJ_EP,type=['pca','tsne','umap','qc'])
		expand(path.join(config['dir']['plot'],'sctNorm','{sample}.sctNorm.clust.qc.jpg'),sample=GEJ_EP),
		expand(path.join(config['dir']['plot'],'rmLCsctNorm','{sample}.rmLCsctNorm.clust.qc.jpg'),sample=GEJ_EP)
	output:
		path.join(config['dir']['log'],'norm_plot_GEJ_EP.finish')
	shell:
		"touch {output}"


rule norm_plot_GEJ_IM:
	input:
		#expand(path.join(config['dir']['plot'],'sctNorm','{sample}.sctNorm.clust.{type}.jpg'),sample=GEJ_IM,type=['pca','tsne','umap','qc']),	
		#expand(path.join(config['dir']['plot'],'rmLCsctNorm','{sample}.rmLCsctNorm.clust.{type}.jpg'),sample=GEJ_IM,type=['pca','tsne','umap','qc'])
		expand(path.join(config['dir']['plot'],'sctNorm','{sample}.sctNorm.clust.qc.jpg'),sample=GEJ_IM),
		expand(path.join(config['dir']['plot'],'rmLCsctNorm','{sample}.rmLCsctNorm.clust.qc.jpg'),sample=GEJ_IM)
	output:
		path.join(config['dir']['log'],'norm_plot_GEJ_IM.finish')
	shell:
		"touch {output}"


rule norm_plot_ESCC_EP:
	input:
		#expand(path.join(config['dir']['plot'],'sctNorm','{sample}.sctNorm.clust.{type}.jpg'),sample=ESCC_EP,type=['pca','tsne','umap','qc']),
		#expand(path.join(config['dir']['plot'],'rmLCsctNorm','{sample}.rmLCsctNorm.clust.{type}.jpg'),sample=ESCC_EP,type=['pca','tsne','umap','qc'])
		expand(path.join(config['dir']['plot'],'sctNorm','{sample}.sctNorm.clust.qc.jpg'),sample=ESCC_EP),
		expand(path.join(config['dir']['plot'],'rmLCsctNorm','{sample}.rmLCsctNorm.clust.qc.jpg'),sample=ESCC_EP)
	output:
		path.join(config['dir']['log'],'norm_plot_ESCC_EP.finish')
	shell:
		"touch {output}"


rule norm_plot_ESCC_IM:
	input:
		#expand(path.join(config['dir']['plot'],'sctNorm','{sample}.sctNorm.clust.{type}.jpg'),sample=ESCC_IM,type=['pca','tsne','umap','qc']),
		#expand(path.join(config['dir']['plot'],'rmLCsctNorm','{sample}.rmLCsctNorm.clust.{type}.jpg'),sample=ESCC_IM,type=['pca','tsne','umap','qc'])
		expand(path.join(config['dir']['plot'],'sctNorm','{sample}.sctNorm.clust.qc.jpg'),sample=ESCC_IM),
		expand(path.join(config['dir']['plot'],'rmLCsctNorm','{sample}.rmLCsctNorm.clust.qc.jpg'),sample=ESCC_IM)
	output:
		path.join(config['dir']['log'],'norm_plot_ESCC_IM.finish')
	shell:
		"touch {output}"


rule plot_all:
	input:
		expand(path.join(config['dir']['log'],'norm_plot_{cancer}_{celltype}.finish'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])

