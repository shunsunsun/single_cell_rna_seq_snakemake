NUM_OF_SUBSAMPLE = config['stable']['num_of_subsample']
ks = config['stable']['subsample_k'].strip().split(',')
#ks = config['stable']['subsample_k']
pcs = config['stable']['subsample_pc'].strip().split(',')
#pcs = config['stable']['subsample_pc']
resoluts = config['stable']['subsample_resolution'].strip()
#resoluts = config['stable']['subsample_resolution']
input_prefix = config['stable']['input_seurat']
cs=['Monocytes','TCells','BCells','ImmuneCells','StromalCells','EpithelialCells']
#cs=['EpithelialCells']

SUBSAMPLE_K_PC = expand(path.join(config['dir']['data'],input_prefix+'_clustStab','subsample',"subsample_k_{k}_PC_{pc}_round_{run_id}.rds"), \
	k = ks, pc = pcs, run_id = range(NUM_OF_SUBSAMPLE))

FULLSAMPLE_K_PC = expand(path.join(config['dir']['data'],input_prefix+'_clustStab','full_sample_preprocess',"full_sample_k_{k}_PC_{pc}.rds"), \
	k = ks, pc = pcs)

GATHER_FULL = expand(path.join(config['dir']['data'],input_prefix+'_clustStab','gather_full_sample.rds'))  
GATHER_SUB = expand(path.join(config['dir']['data'],input_prefix+'_clustStab','gather_subsample.rds'))

#TARGETS = []
#TARGETS.extend(SUBSAMPLE_K_PC)
#TARGETS.append(GATHER_FULL)
#TARGETS.append(GATHER_SUB)

#rule clust_stab_eval:
#	input: TARGETS
#	output: path.join(config['dir']['log'],"clust_stab_sampling.finish")
#	shell:
#		"touch {output}"

combined_method="sctNorm_BatchHmy"
#combined_method="rc_glmpca"

## the full data set, preprocessing using a set of k, resolution and PC
rule full_sample_preprocess:
	#input: path.join(config['dir']['data'], input_prefix + ".rds")
	#input: path.join(config['dir']['data'], input_prefix + '_clustStab', '{celltype}_sctNorm_BatchHmy.rds')
	input: path.join(config['dir']['data'], 'Combined', '{celltype}_'+combined_method+'.rds')
	#output: path.join(config['dir']['data'],input_prefix+'_clustStab','full_sample_preprocess', "full_sample_k_{k}_PC_{pc}.rds")
	#output: path.join(config['dir']['data'],input_prefix+'_clustStab','{celltype}_full_sample_preprocess','full_sample_k_{k}_PC_{pc}.rds')
	output: path.join(config['dir']['data'], 'Combined', '{celltype}_'+combined_method+'_full_sample_preprocess', 'full_sample_k_{k}_PC_{pc}.rds')
	#log: path.join(config['dir']['log'], input_prefix+"_full_sample_k_{k}_PC_{pc}.log")
	params: res=resoluts,
		reduc=config['stable']['stability_dimReduct'],
		PreprocessSubsetData_pars=config['stable']['PreprocessSubsetData_subsample_pars']
	threads: 36
	#message: "preprocessing original full seurat object using k of {wildcards.k} and {wildcards.pc} PCs with {threads} threads"
	script: "../scripts/clustStab/preprocess.R"


rule gather_full_sample_preprocess:
	#input: rds = FULLSAMPLE_K_PC
	#input: rds = lambda wildcard: expand(path.join(config['dir']['data'],input_prefix+'_clustStab','{{celltype}}_full_sample_preprocess','full_sample_k_{k}_PC_{pc}.rds'), k=ks, pc = pcs)
	input: rds = lambda wildcard: expand(path.join(config['dir']['data'],'Combined', '{{celltype}}_'+combined_method+'_full_sample_preprocess', 'full_sample_k_{k}_PC_{pc}.rds'),k=ks,pc=pcs)
	#output: GATHER_FULL
	#output: path.join(config['dir']['data'],input_prefix+'_clustStab','{celltype}_full_sample_preprocess','gather_full_sample.rds')
	output: path.join(config['dir']['data'],'Combined', '{celltype}_'+combined_method+'_full_sample_preprocess','gather_full_sample.rds')
	params: res=resoluts
	#log: path.join(config['dir']['log'], input_prefix+"_full_sample_gather_idents.log")
	threads: 36
	#message: "gathering full sample clustering results"
	script: "../scripts/clustStab/gather_fullsample.R"


## subsample e.g. 80% of the cells and re-do the clustering for n times
rule subsample_cluster:
	#input: path.join(config['dir']['data'],input_prefix+'_clustStab','full_sample_preprocess', "full_sample_k_{k}_PC_{pc}.rds")
	#input: path.join(config['dir']['data'],input_prefix+'_clustStab','{celltype}_full_sample_preprocess','full_sample_k_{k}_PC_{pc}.rds')
	input: path.join(config['dir']['data'], 'Combined', '{celltype}_'+combined_method+'_full_sample_preprocess', 'full_sample_k_{k}_PC_{pc}.rds')
	#output: path.join(config['dir']['data'],input_prefix+'_clustStab','subsample',"subsample_k_{k}_PC_{pc}_round_{run_id}.rds")
	#output: path.join(config['dir']['data'],input_prefix+'_clustStab','{celltype}_subsample', 'subsample_k_{k}_PC_{pc}_round_{run_id}.rds')
	output: path.join(config['dir']['data'],'Combined', '{celltype}_'+combined_method+'_subsample', 'subsample_k_{k}_PC_{pc}_round_{run_id}.rds')
	#log: path.join(config['dir']['log'],input_prefix+"subsample_k_{k}_PC_{pc}_round_{run_id}.log")
	params: rate = config['stable']['subsample_rate'],
		regressNum=config['norm']['regressNum'],
		res = resoluts,
		reduc = config['stable']['stability_dimReduct'],
		PreprocessSubsetData_pars = config['stable']['PreprocessSubsetData_subsample_pars']
	threads: 36
	#message: "subsampling {params.rate} from the full data set, recluster using k of {wildcards.k} and {wildcards.pc} PCs for round {wildcards.run_id} using {threads} threads"
	script: "../scripts/clustStab/subsample.R"


## gather the subsampled and reclustered cell idents
rule gather_subsample:
	#input: rds = SUBSAMPLE_K_PC
	#input: rds = lambda wildcards: expand(path.join(config['dir']['data'],input_prefix+'_clustStab','{{celltype}}_subsample','subsample_k_{k}_PC_{pc}_round_{run_id}.rds'), k = ks, pc = pcs, run_id = range(NUM_OF_SUBSAMPLE))
	input: rds = lambda wildcards: expand(path.join(config['dir']['data'],'Combined', '{{celltype}}_'+combined_method+'_subsample', 'subsample_k_{k}_PC_{pc}_round_{run_id}.rds'),k = ks, pc = pcs, run_id = range(NUM_OF_SUBSAMPLE))
	#output: GATHER_SUB
	#output: path.join(config['dir']['data'],input_prefix+'_clustStab','{celltype}_subsample','gather_subsample.rds')
	output: path.join(config['dir']['data'],'Combined','{celltype}_'+combined_method+'_subsample','gather_subsample.rds')
	params: res=resoluts
	#log: path.join(config['dir']['log'],input_prefix+"_gather_subsample.log")
	threads: 36
	#message: "gathering clustering results for subsamples"
	script: "../scripts/clustStab/gather_subsample.R"


rule summarize_res:
	input:
		#full = GATHER_FULL,
		#sub = GATHER_SUB
		#full = path.join(config['dir']['data'],input_prefix+'_clustStab','{celltype}_full_sample_preprocess','gather_full_sample.rds'),
		#sub = path.join(config['dir']['data'],input_prefix+'_clustStab','{celltype}_subsample','gather_subsample.rds')
		full=path.join(config['dir']['data'],'Combined', '{celltype}_'+combined_method+'_full_sample_preprocess','gather_full_sample.rds'),
		sub=path.join(config['dir']['data'],'Combined','{celltype}_'+combined_method+'_subsample','gather_subsample.rds')
	output:
		#path.join(config['dir']['plot'],'clustStab',input_prefix+'_scatterPlot1.pdf'),
		#path.join(config['dir']['plot'],'clustStab',input_prefix+'_{celltype}_scatterPlot1_v2.pdf'),
		path.join(config['dir']['plot'],'clustStab',combined_method+'_{celltype}_scatterPlot1.pdf'),
		#path.join(config['dir']['plot'],'clustStab',input_prefix+'_scatterPlot3.pdf')
		#path.join(config['dir']['plot'],'clustStab',input_prefix+'_{celltype}_scatterPlot3_v2.pdf')
		path.join(config['dir']['plot'],'clustStab',combined_method+'_{celltype}_scatterPlot3.pdf')
	threads: 36
	script: "../scripts/clustStab/analyze_plot.R"


rule summarize_all:
	input:
		#expand(path.join(config['dir']['plot'],'clustStab',input_prefix+'_{celltype}_scatterPlot1_v2.pdf'), celltype=cs),
		#expand(path.join(config['dir']['plot'],'clustStab',input_prefix+'_{celltype}_scatterPlot3_v2.pdf'), celltype=cs)
		expand(path.join(config['dir']['plot'],'clustStab',combined_method+'_{celltype}_scatterPlot1.pdf'), celltype=cs),
		expand(path.join(config['dir']['plot'],'clustStab',combined_method+'_{celltype}_scatterPlot3.pdf'), celltype=cs)
	output:
		path.join(config['dir']['log'], 'summarize_all.finish')
	shell:
		"touch {output}"
