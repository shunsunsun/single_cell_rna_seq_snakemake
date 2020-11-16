reduc={'GEJ':'glmpca10','ESCC':'glmpca20'}
SingleRref={'EP':'hpca_be','IM':'hpca_be_dice_nh_mi'}

rule ident_by_SingleR:
        input:
                path.join(config['dir']['data'],'{cancer}_{celltype}_pipeCompflted_clustered.rds')
        output:
                path.join(config['dir']['data'],'{cancer}_{celltype}_pipeCompflted_clustered_SingleRann.rds'),
		path.join(config['dir']['plot'],'cellident','{cancer}_{celltype}_SingleRlabelTSNEUMAP.jpg')
	params:
		interimRes=path.join(config['dir']['data'],'{cancer}_{celltype}_SingleRpred.rda'),
		label=config['annot']['SingleRlabel'],
		ref=lambda wildcards: SingleRref[wildcards.celltype],
		refdir=path.join(config['dir']['resources'],'SingleRdata'),
		reducPrefix=lambda wildcards: reduc[wildcards.cancer]
	threads:
		40
	script:
                "../scripts/step10_cellident_byRef.R"


rule ident_by_scHCL:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_pipeCompflted_clustered_SingleRann.rds')
	output:
		path.join(config['dir']['data'],'{cancer}_{celltype}_SingleR_scHCL_compare.rds')
	params:
		interimRes=path.join(config['dir']['data'],'{cancer}_{celltype}_scHCLres.rda')
	threads:
		3
	script:
		"../scripts/step10_cellident_byscHCL.R"


rule ident_by_scHCL_SingleR:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_pipeCompflted_clustered_SingleRann.rds'),
		path.join(config['dir']['data'],'{cancer}_{celltype}_SingleR_scHCL_compare.rds')
	output:
		path.join(config['dir']['plot'],'cellident','{cancer}_{celltype}_scHCL_sankeyPlot2.jpg')
	params:
		interimRes=path.join(config['dir']['data'],'{cancer}_{celltype}_scHCL_SingleRpred.rda'),
                ref=path.join(config['dir']['resources'],'SingleRdata','scHCL_ref.rds'),
                reducPrefix=lambda wildcards: reduc[wildcards.cancer]
	threads:
		40
	script:
		"../scripts/step10_cellident_bySCref.R"


rule ident_by_markers:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_pipeCompflted_clustered_SingleRann.rds')
	output:
		path.join(config['dir']['data'],'{cancer}_{celltype}_MarkerPred.txt')
	params:
		outPlot=path.join(config['dir']['plot'],'cellident'),
		reducPrefix=lambda wildcards: reduc[wildcards.cancer]
	script:
		"../scripts/step10_cellident_byMarkers.R"


def get_input(wildcards):
	wanted = []
	if config['annot']['SingleR']==True:
		wanted.extend(expand(path.join(config['dir']['plot'],'cellident','{cancer}_{celltype}_SingleRlabelTSNEUMAP.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['annot']['scHCL']==True:
		wanted.extend(expand(path.join(config['dir']['plot'],'cellident','{cancer}_{celltype}_scHCL_sankeyPlot2.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['annot']['Marker']==True:
		wanted.extend(expand(path.join(config['dir']['data'],'{cancer}_{celltype}_MarkerPred.txt'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	return wanted


rule ident_all:
	input:
		get_input
	output:
		path.join(config['dir']['log'],'cell_ident.finish')
	shell:
		"touch {output}"
