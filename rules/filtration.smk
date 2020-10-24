rule filter_hardThresh:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_merged.rds')
	params:
		outdir=config['dir']['data'],
		nUMI_lo=config['qc']['nUMI_lo'],
		nUMI_up=config['qc']['nUMI_up'],
		nGene_lo=config['qc']['nGene_lo'],
		nGene_up=config['qc']['nGene_up'],
		log10GenesPerUMI=config['qc']['log10GenesPerUMI'],
		mitoCountRatio=config['qc']['mitoCountRatio'],
		nCell=config['qc']['nCell']
	output:
        	path.join(config['dir']['data'],'{cancer}_{celltype}_hardflted.rds')
	script:
		"../scripts/filter_hardThresh.R"


rule filter_outlier:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_merged.rds')
	params:
		all=config['qc']['all_sample'],
		lowqual_list=config['manifest']['lowqual_sample'],
		plotDir=path.join(config['dir']['qc'],'outlyflted'),
		outdir=config['dir']['data'],
		nCell=config['qc']['nCell']
	output:
		rds=path.join(config['dir']['data'],'{cancer}_{celltype}_outlyflted.rds'),
		stats=path.join(config['dir']['log'],'{cancer}_{celltype}_outlyflted_stat.tsv')
	script:
		"../scripts/step5_removeOutliers.R"


rule filter_pipeComp:
	input:
		path.join(config['dir']['data'],'{cancer}_{celltype}_rmDbl_sce.rds')
	params:
		all=config['qc']['all_sample'],
		lowqual_list=config['manifest']['lowqual_sample'],
		outdir=config['dir']['data'],
		nCell=config['qc']['nCell'],
		mtRatio=config['qc']['mitoCountRatio']
	output:
		rds=path.join(config['dir']['data'],'{cancer}_{celltype}_pipeCompflted.rds'),
	script:
		"../scripts/step5_removeOutliers_pipeComp.R"



if 'hard' in config['qc']['strategy']:
	label='hardflted'
if 'outly' in config['qc']['strategy']:
	label='outlyflted'
if 'pipe' in config['qc']['strategy']:
	label='pipeCompflted'


rule filter_all:
	input:
		expand(path.join(config['dir']['data'],'{cancer}_{celltype}_'+label+'.rds'),cancer=['ESCC','GEJ'],celltype=['EP','IM'])
	output:
		path.join(config['dir']['log'],label + '_filter.finished')
	shell:
		"touch {output}"
