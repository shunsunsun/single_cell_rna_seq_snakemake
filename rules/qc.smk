rule qc_plot:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_merged.rds')
	params:
		config['dir']['qc']
	output:
		path.join(config['dir']['qc'],'{cancer}_{celltype}_VlnPlot_mito_ribo_perc.jpg')
		path.join(config['dir']['qc'],'{cancer}_{celltype}_VlnPlot_nCount_5Kmax.jpg')
		path.join(config['dir']['qc'],'{cancer}_{celltype}_VlnPlot_nFeature.jpg')
		path.join(config['dir']['qc'],'{cancer}_{celltype}_VlnPlot_nGenePerUMI.jpg')
		path.join(config['dir']['qc'],'{cancer}_{celltype}_ScattPlot_nCount_mitoRatio.jpg')
		path.join(config['dir']['qc'],'{cancer}_{celltype}_ScattPlot_nCount_riboRatio.jpg')
		path.join(config['dir']['qc'],'{cancer}_{celltype}_ScattPlot_nCount_nFeature.jpg')
		path.join(config['dir']['qc'],'{cancer}_{celltype}_ScattPlot_nCount_nFeature_perSample.jpg')
	script:
		"../scripts/qc_plot.R"


rule qc_plot_all:
	input:
		expand(path.join(config['dir']['qc'],'{cancer}_{celltype}_VlnPlot_mito_ribo_perc.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
		expand(path.join(config['dir']['qc'],'{cancer}_{celltype}_VlnPlot_nCount_5Kmax.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
		expand(path.join(config['dir']['qc'],'{cancer}_{celltype}_VlnPlot_nFeature.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
		expand(path.join(config['dir']['qc'],'{cancer}_{celltype}_VlnPlot_nGenePerUMI.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
		expand(path.join(config['dir']['qc'],'{cancer}_{celltype}_ScattPlot_nCount_mitoRatio.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
		expand(path.join(config['dir']['qc'],'{cancer}_{celltype}_ScattPlot_nCount_riboRatio.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
		expand(path.join(config['dir']['qc'],'{cancer}_{celltype}_ScattPlot_nCount_nFeature.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
		expand(path.join(config['dir']['qc'],'{cancer}_{celltype}_ScattPlot_nCount_nFeature_perSample.jpg'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['log'],"qc_plot.finish")
	shell:
		"touch {output}"


rule qc_discard:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_merged.rds')
	params:
		all=config['qc']['all_sample']
		lowqual_list=config['manifest']['lowqual_sample']
		outlying=config['qc']['cal_outlyingness']
		outdir=config['dir']['qc']
	output:
		path.join(config['dir']['qc'],'{cancer}_{celltype}_scater_outlier.jpg')
	script:
		"../scripts/qc_discard.R"



rule qc_discard_all:
	input:
		
