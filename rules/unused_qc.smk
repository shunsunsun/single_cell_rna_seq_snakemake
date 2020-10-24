rule qc_metrics:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_merged_noQCmetrics.rds')
	output:
		path.join(config['dir']['data'],'{cancer}_{celltype}_merged.rds')
	threads:
		4
	params:
		house=path.join(config['dir']['resources'],'HRT_atlas_Human_HousekeepingGenes.RData')		
	script:
		"../scripts/step3_calculateQCmetrics.R"


rule qc_metrics_all:
	input:
		expand(path.join(config['dir']['data'],'{cancer}_{celltype}_merged.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],'calculate_metrics.finish')
	shell:
		"touch {output}"



label=config['qc']['target']
## return a list of wanted plots based on the config file
def wanted_qc_plot():
	wanted = []
	if config['plot']['bar_nCell'] == true:
		wanted.extend(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_BarPlot_nCellperSample.jpg'))
	if config['plot']['vln_mito'] == true:
		wanted.extend(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_mito.jpg'))
	if config['plot']['vln_ribo'] == true:
		wanted.extend(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_ribo.jpg'))
	if config['plot']['vln_hk'] == true:
		wanted.extend(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_hk.jpg'))
	if config['plot']['vln_hb'] == true:
		wanted.extend(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_hb.jpg'))
	if config['plot']['vln_nCount'] == true:
		wanted.extend(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nCount_5Kmax.jpg'))
	if config['plot']['vln_nFeature'] == true:
		wanted.extend(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nFeature.jpg'))
	if config['plot']['vln_complexity'] == true:
		wanted.extend(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nGenePerUMI.jpg'))
	if config['plot']['scatter_nCount_mito'] == true:
		wanted.extend(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_mito.jpg'))
	if config['plot']['scatter_nCount_nFeature_colSample'] == true:
		wanted.extend(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature.jpg'))
	if config['plot']['scatter_nCount_nFeature_colHb'] == true:
		wanted.extend(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_colHb.jpg'))
	if config['plot']['scatter_nCount_nFeature_perSample'] == true:
		wanted.extend(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_perSample.jpg'))
	return wanted


rule qc_plot:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_' + label + '.rds')
	params:
		outdir=path.join(config['dir']['qc'],label),
		nCell=config['qc']['nCell'],
		bar_nCell=config['plot']['bar_nCell'],
		vln_mito=config['plot']['vln_mito'],
		vln_ribo=config['plot']['vln_ribo'],
		vln_hk=config['plot']['vln_hk'],
		vln_hb=config['plot']['vln_hb'],
		vln_nCount=config['plot']['vln_nCount'],
		vln_nFeature=config['plot']['vln_nFeature'],
		vln_complexity=config['plot']['vln_complexity'],
		scatter_nCount_mito=config['plot']['scatter_nCount_mito'],
		scatter_nCount_nFeature_colSample=config['plot']['scatter_nCount_nFeature_colSample'],
		scatter_nCount_nFeature_colHb=config['plot']['scatter_nCount_nFeature_colHb'],
		scatter_nCount_nFeature_perSample=config['plot']['scatter_nCount_nFeature_perSample']
	output:
		wanted_qc_plot
	script:
		"../scripts/step4_QCmetricsPlot.R"


def get_input(wildcards):
	wanted = []
	if config['plot']['bar_nCell'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_BarPlot_nCellperSample.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_mito'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_mito.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_ribo'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_ribo.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_hk'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_hk.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_hb'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_hb.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_nCount'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nCount_5Kmax.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_nFeature'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nFeature.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_complexity'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nGenePerUMI.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['scatter_nCount_mito'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_mito.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['scatter_nCount_nFeature_colSample'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['scatter_nCount_nFeature_colSample'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_colSample.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['scatter_nCount_nFeature_colHb'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_colHb.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['scatter_nCount_nFeature_perSample'] == true:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_perSample.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	return wanted


rule qc_plot_all:
	input:
		get_input
	output:
		path.join(config['dir']['log'],'qc_plot_'+label+'.finish')
	shell:
		"touch {output}"


if config['qc']['all_sample']:
	label2='_discard_allbased'  ##based on all samples
else:
	label2='_discard_hqbased'   ##based on high-quality samples


rule qc_discard:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_' + label + '.rds')
	params:
		all=config['qc']['all_sample'],
		lowqual_list=config['manifest']['lowqual_sample'],
		outlying=config['qc']['cal_outlyingness'],
		outdir=path.join(config['dir']['qc'],label),
	output:
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_nCount.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_nFeature.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_mitoPerc.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_nCount_mitoPerc.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_nCount_mitoPerc_perSample.jpg')
	script:
		"../scripts/qc_discard_diagnoPlot.R"


rule qc_discard_all:
	input:
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2+ '_nCount.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_nFeature.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_mitoPerc.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_nCount_mitoPerc.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_nCount_mitoPerc_perSample.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],'qc_' + label + label2 + '.finish')
	shell:
		"touch {output}"
		
