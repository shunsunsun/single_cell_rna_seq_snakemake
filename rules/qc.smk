label=config['qc']['target']

rule qc_plot:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_' + label + '.rds')
	params:
		outdir=path.join(config['dir']['qc'],label),
		house=path.join(config['dir']['resources'],'HRT_atlas_Human_HousekeepingGenes.RData'),
		nCell=config['qc']['nCell']
	output:
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}_BarPlot_nCellperSample.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_ribo_hk_perc.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_mito_perc.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nCount_5Kmax.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nFeature.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nGenePerUMI.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_mitoRatio.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_riboRatio.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature.jpg'),
		path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_perSample.jpg')
	script:
		"../scripts/qc_plot.R"


rule qc_plot_all:
	input:
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_BarPlot_nCellperSample.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_ribo_hk_perc.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_mito_perc.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nCount_5Kmax.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nFeature.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_VlnPlot_nGenePerUMI.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_mitoRatio.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_riboRatio.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM']),
		expand(path.join(config['dir']['qc'],label,'{cancer}_{celltype}_ScattPlot_nCount_nFeature_perSample.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
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
		
