rule createSeuratObject:
	input:
		path.join(config['dir']['rawdata'],'{sample}','GRCh38','genes.tsv'),
		path.join(config['dir']['rawdata'],'{sample}','GRCh38','barcodes.tsv'),
		path.join(config['dir']['rawdata'],'{sample}','GRCh38','matrix.mtx')
	params:
		dir=path.join(config['dir']['rawdata'],'{sample}','GRCh38'),
		id='{sample}'
	output:
		path.join(config['dir']['data'],'initSeuObj','{sample}.initSeu.rds')
	script:
		"../scripts/step1_sample2SeuratObj.R"		


rule merge_GEJ_EP:
	input:
		objs=expand(path.join(config['dir']['data'],'initSeuObj','{sample}.initSeu.rds'),sample=GEJ_EP),
		list=path.join('config','samples_GEJ_EP.txt')
	output:
		path.join(config['dir']['data'],'GEJ_EP_merged_noQCmetrics.rds')
	params:
		indir=config['dir']['data'],
		projectName="gej_ep",
		suffix=".seuratobj.rds"
	script:
		"../scripts/step2_mergeSamples.R"


rule merge_GEJ_IM:
	input:
		objs=expand(path.join(config['dir']['data'],'initSeuObj','{sample}.initSeu.rds'),sample=GEJ_IM),
		list=path.join('config','samples_GEJ_IM.txt')
	output:
		path.join(config['dir']['data'],'GEJ_IM_merged_noQCmetrics.rds')
	params:
		indir=config['dir']['data'],
		projectName="gej_im",
		suffix=".seuratobj.rds"
	script:
		"../scripts/step2_mergeSamples.R"


rule merge_ESCC_EP:
	input:
		objs=expand(path.join(config['dir']['data'],'initSeuObj','{sample}.initSeu.rds'),sample=ESCC_EP),
		list=path.join('config','samples_ESCC_EP.txt')
	output:
		path.join(config['dir']['data'],'ESCC_EP_merged_noQCmetrics.rds')
	params:
		indir=config['dir']['data'],
		projectName="escc_ep",
		suffix=".seuratobj.rds"
	script:
		"../scripts/step2_mergeSamples.R"


rule merge_ESCC_IM:
	input:
		objs=expand(path.join(config['dir']['data'],'initSeuObj','{sample}.initSeu.rds'),sample=ESCC_IM),
		list=path.join('config','samples_ESCC_IM.txt')
	output:
		path.join(config['dir']['data'],'ESCC_IM_merged_noQCmetrics.rds')
	params:
		indir=config['dir']['data'],
		projectName="escc_im",
		suffix=".seuratobj.rds"
	script:
		"../scripts/step2_mergeSamples.R"
	

rule merge_all:
	input:
		expand(path.join(config['dir']['data'],'{cancer}_{celltype}_merged.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],'merge_samples.finish')
	shell:
		"touch {output}"


rule merge_ESCC:
	input:
		path.join('config','merge_escc_samples.txt')
	output:
		path.join(config['dir']['data'],'ESCC_QCed.rds')
	params:
		"ESCC"
	script:
		"../scripts/step2_mergeDataSets.R"


rule merge_GEJ:
	input:
		path.join('config','merge_gej_samples.txt')
	output:
		path.join(config['dir']['data'],'GEJ_QCed.rds')
	params:
		"GEJ"
	script:
		"../scripts/step2_mergeDataSets.R"


rule merge_IM:
	input:
		path.join('config','merge_immune_samples.txt')
	output:
		path.join(config['dir']['data'],'IM_QCed.rds')
	params:
		"ESCC_GEJ_IM"
	script:
		"../scripts/step2_mergeDataSets.R"


rule merge_datasets_all:
	input:
		expand(path.join(config['dir']['data'], '{dataset}_QCed.rds'),dataset=['ESCC','GEJ','IM'])
	output:
		path.join(config['dir']['log'],'merge_datasets.finish')
	shell:
		"touch {output}"
