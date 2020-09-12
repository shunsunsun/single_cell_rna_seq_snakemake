rule createSeuratObject:
	input:
		path.join(config['dir']['rawdata'],'{sample}','GRCh38','genes.tsv'),
		path.join(config['dir']['rawdata'],'{sample}','GRCh38','barcodes.tsv'),
		path.join(config['dir']['rawdata'],'{sample}','GRCh38','matrix.mtx')
	params:
		dir=path.join(config['dir']['rawdata'],'{sample}','GRCh38'),
		id='{sample}'
	output:
		path.join(config['dir']['data'],'{sample}.seuratobj.rds')
	script:
		"../scripts/create_seurat_obj.R"		


rule merge_GEJ_EP:
	input:
		objs=expand(path.join(config['dir']['data'],'{sample}.seuratobj.rds'),sample=GEJ_EP),
		list=path.join('config','samples_GEJ_epis.txt')
	output:
		path.join(config['dir']['data'],'GEJ_EP_merged.rds')
	params:
		config['dir']['data']
	script:
		"../scripts/merge_cells.R"


rule merge_GEJ_IM:
	input:
		objs=expand(path.join(config['dir']['data'],'{sample}.seuratobj.rds'),sample=GEJ_IM),
		list=path.join('config','samples_GEJ_immun.txt')
	output:
		path.join(config['dir']['data'],'GEJ_IM_merged.rds')
	params:
		config['dir']['data']
	script:
		"../scripts/merge_cells.R"


rule merge_ESCC_EP:
	input:
		objs=expand(path.join(config['dir']['data'],'{sample}.seuratobj.rds'),sample=ESCC_EP),
		list=path.join('config','samples_ESCC_epis.txt')
	output:
		path.join(config['dir']['data'],'ESCC_EP_merged.rds')
	params:
		config['dir']['data']
	script:
		"../scripts/merge_cells.R"


rule merge_ESCC_IM:
	input:
		objs=expand(path.join(config['dir']['data'],'{sample}.seuratobj.rds'),sample=ESCC_IM),
		list=path.join('config','samples_ESCC_immun.txt')
	output:
		path.join(config['dir']['data'],'ESCC_IM_merged.rds')
	params:
		config['dir']['data']
	script:
		"../scripts/merge_cells.R"
	

rule merge_all:
	input:
		expand(path.join(config['dir']['data'],'{cancer}_{celltype}_merged.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],'merge_cell.finish')
	shell:
		"touch {output}"
