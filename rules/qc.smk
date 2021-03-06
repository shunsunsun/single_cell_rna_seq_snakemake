rule qc_metrics:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_merged_noQCmetrics.rds')
	output:
		path.join(config['dir']['data'],'{cancer}_{celltype}_merged.rds')
	threads: 4
	params:
		path.join(config['dir']['resources'],'HRT_atlas_Human_HousekeepingGenes.RData')
	script:
		"../scripts/step3_calculateQCmetrics.R"


rule qc_metrics_scater:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_merged_noQCmetrics.rds')
	output:
		path.join(config['dir']['data'],'{cancer}_{celltype}_rmDbl_sce.rds')
	threads: 4
	script:
		"../scripts/step3_calculateQCmetrics_scater.R"


rule qc_metrics_all:
	input:
		expand(path.join(config['dir']['data'],'{cancer}_{celltype}_rmDbl_sce.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],'scater_qc_metrics.finish')
	shell:
		"touch {output}"



label=config['qc']['target']


rule qc_plot:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_' + label + '.rds')
	params:
		outdir=path.join(config['dir']['qc'],label),
		#nCell=config['qc']['nCell'],
                #bar_nCell=path.join(config['plot'],'bar_nCell'),
                #vln_mito=path.join(config['plot'],'vln_mito'),
                #vln_ribo=path.join(config['plot'],'vln_ribo'),
                #vln_hk=path.join(config['plot'],'vln_hk'),
                #vln_hb=path.join(config['plot'],'vln_hb'),
                #vln_nCount=path.join(config['plot'],'vln_nCount'),
                #vln_nFeature=path.join(config['plot'],'vln_nFeature'),
                #vln_complexity=path.join(config['plot'],'vln_complexity'),
                #scatter_nCount_mito=path.join(config['plot'],'scatter_nCount_mito'),
                #scatter_nCount_nFeature_colSample=path.join(config['plot'],'scatter_nCount_nFeature_colSample'),
                #scatter_nCount_nFeature_colHb=path.join(config['plot'],'scatter_nCount_nFeature_colHb'),
                #scatter_nCount_nFeature_perSample=path.join(config['plot'],'scatter_nCount_nFeature_perSample')
	output:
		#path.join(config['dir']['qc'],label,'{cancer}_{celltype}_BarPlot_nCellperSample.jpg'),
		#path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_ribo.jpg'),
		#path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_mito.jpg'),
		#path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_hk.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_hb.jpg'),
		#path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nCount_5Kmax.jpg'),
		#path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nFeature.jpg'),
		#path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nGenePerUMI.jpg'),
		#path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_mito.jpg'),
		#path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature.jpg'),
		#path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_colSample.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_colHb.jpg')
	script:
		"../scripts/step4_plot_QCmetrics.R"


def get_qc_input(wildcards):
	wanted = []
	if config['plot']['bar_nCell'] == True:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_BarPlot_nCellperSample.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_mito'] == True:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_mito.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_ribo'] == True:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_ribo.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_hk'] == True:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_hk.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_hb'] == True:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_hb.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_nCount'] == True:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nCount_5Kmax.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_nFeature'] == True:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nFeature.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['vln_complexity'] == True:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nGenePerUMI.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['scatter_nCount_mito'] == True:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_mito.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['scatter_nCount_nFeature_colSample'] == True:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['scatter_nCount_nFeature_colHb'] == True:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_colHb.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	if config['plot']['scatter_nCount_nFeature_perSample'] == True:
		wanted.extend(expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_perSample.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']))
	return wanted


rule qc_plot_all:
	input:
		get_qc_input
		#expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_BarPlot_nCellperSample.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		#expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_ribo.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		#expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_mito.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		#expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_hk.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		#expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nCount_5Kmax.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		#expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nFeature.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		#expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nGenePerUMI.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		#expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_mito.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		#expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		#expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_perSample.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		#expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_colHb.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		#expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_hb.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],'qc_plot_'+label+'.finish')
	shell:
		"touch {output}"


if config['qc']['all_sample']:
	label2='_discard_allbased'  ##based on all samples
else:
	label2='_discard_hqbased'   ##based on high-quality samples


rule plot_discard:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_' + label + '.rds')
	params:
		all=config['qc']['all_sample'],
		lowqual_list=config['manifest']['lowqual_sample'],
		#outlying=config['qc']['cal_outlyingness'],
		outdir=path.join(config['dir']['qc'],label),
	output:
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_nCount.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_nFeature.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_mitoPerc.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_nCount_mitoPerc.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}' + label2 + '_nCount_mitoPerc_perSample.jpg')
	script:
		"../scripts/step5b_plot_discard.R"


rule plot_discard_all:
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
		
