if 'hard' in config['qc']['strategy']:
	inData='hardflted'
if 'outly' in config['qc']['strategy']:
	inData='outlyflted'
if 'pipe' in config['qc']['strategy']:
	inData='pipeCompflted'


rule seurat3_integration:
	input:
		list=path.join('config','samples_{cancer}_{celltype}.txt')
	output:
		rds=path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_IntegBySeurat3.rds'),
		#pca=path.join(config['dir']['plot'],'integrate','{cancer}_{celltype}_'+inData+'_IntegBySeurat3_pca.jpg'),
		#tsne=path.join(config['dir']['plot'],'integrate','{cancer}_{celltype}_'+inData+'_IntegBySeurat3_tsne.jpg'),
		#umap=path.join(config['dir']['plot'],'integrate','{cancer}_{celltype}_'+inData+'_IntegBySeurat3_umap.jpg'),
	params:
		indir=path.join(config['dir']['data'],'sctNormSeuObj_'+inData),
		infile="sctNormSeu_" + inData,
		nPC=config['integrate']['seurat3_nPC'],
		nFeature=config['integrate']['seurat3_nFeature'],
		anchor=config['integrate']['seurat3_anchor']		
	threads:
		20
	script:
		"../scripts/step7_integrateSamples_seurat3.R"


def get_integrate_method(wildcards):
	wanted = []
	if config['integrate']['runSeurat3'] == True:
		wanted.extend(expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+inData+'_IntegBySeurat3.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	return wanted


rule integrate_all:
	input:
		get_integrate_method
	output:
		path.join(config['dir']['log'],'integrate_all.finish')
	shell:
		"touch {output}"


