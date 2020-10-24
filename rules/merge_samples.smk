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
		config['dir']['data']
	script:
		"../scripts/step2_mergeSamples.R"


rule merge_GEJ_IM:
	input:
		objs=expand(path.join(config['dir']['data'],'initSeuObj','{sample}.initSeu.rds'),sample=GEJ_IM),
		list=path.join('config','samples_GEJ_IM.txt')
	output:
		path.join(config['dir']['data'],'GEJ_IM_merged_noQCmetrics.rds')
	params:
		config['dir']['data']
	script:
		"../scripts/step2_mergeSamples.R"


rule merge_ESCC_EP:
	input:
		objs=expand(path.join(config['dir']['data'],'initSeuObj','{sample}.initSeu.rds'),sample=ESCC_EP),
		list=path.join('config','samples_ESCC_EP.txt')
	output:
		path.join(config['dir']['data'],'ESCC_EP_merged_noQCmetrics.rds')
	params:
		config['dir']['data']
	script:
		"../scripts/step2_mergeSamples.R"


rule merge_ESCC_IM:
	input:
		objs=expand(path.join(config['dir']['data'],'initSeuObj','{sample}.initSeu.rds'),sample=ESCC_IM),
		list=path.join('config','samples_ESCC_IM.txt')
	output:
		path.join(config['dir']['data'],'ESCC_IM_merged_noQCmetrics.rds')
	params:
		config['dir']['data']
	script:
		"../scripts/step2_mergeSamples.R"
	

rule merge_all:
	input:
		expand(path.join(config['dir']['data'],'{cancer}_{celltype}_merged.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],'merge_samples.finish')
	shell:
		"touch {output}"
