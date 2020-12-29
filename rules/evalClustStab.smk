NUM_OF_SUBSAMPLE = config['stable']['num_of_subsample']
ks = config['stable']['subsample_k'].strip().split(',')
pcs = config['stable']['subsample_pc'].strip().split(',')

SUBSAMPLE_K_PC = expand(path.join(config['dir']['data'],'evalClustStab',"subsample/subsample_k_{k}_PC_{pc}_round_{run_id}.rds"), \
	k = ks, pc = pcs, run_id = range(NUM_OF_SUBSAMPLE))

FULLSAMPLE_K_PC = expand(path.join(config['dir']['data'],'evalClustStab',"full_sample_preprocess/full_sample_k_{k}_PC_{pc}.rds"), \
	k = ks, pc = pcs)

TARGETS = []
TARGETS.extend(SUBSAMPLE_K_PC)
TARGETS.append("gather_subsample.rds")
TARGETS.append("gather_full_sample.rds")

rule clust_stab_eval_sampling:
	input: TARGETS
	output: path.join(config['dir']['log'],"clust_stab_sampling.finish")
	shell:
		"touch {output}"


## the full data set, preprocessing using a set of k, resolution and PC
rule full_sample_preprocess:
	input: config['stable']['input_seurat']
	output: path.join(config['dir']['data'],'evalClustStab',"full_sample_preprocess/full_sample_k_{k}_PC_{pc}.rds")
	log: path.join(config['dir']['log'], "full_sample_k_{k}_PC_{pc}.log")
	params: jobname = "full_sample_k_{k}_PC_{pc}",
		resolution = config['stable']['subsample_resolution'].strip(),
		PreprocessSubsetData_pars = config.get("PreprocessSubsetData_subsample_pars", "")
	threads: 12
	message: "preprocessing original full seurat object using k of {wildcards.k} and {wildcards.pc} PCs with {threads} threads"
	script: "../scripts/clustStab/preprocess.R"


rule gather_full_sample_preprocess:
	input: rds = FULLSAMPLE_K_PC
	output: path.join(config['dir']['data'],'evalClustStab',"gather_full_sample.rds")
	log: path.join(config['dir']['log'], "full_sample_gather_idents.log")
	threads: 3
	message: "gathering full sample idents"
	script: "../scripts/clustStab/gather_fullsample.R"


## subsample e.g. 80% of the cells and re-do the clustering for n times
rule subsample_cluster:
	input: path.join(config['dir']['data'],'evalClustStab',"full_sample_preprocess/full_sample_k_{k}_PC_{pc}.rds")
	output: path.join(config['dir']['data'],'evalClustStab',"subsample/subsample_k_{k}_PC_{pc}_round_{run_id}.rds")
	log: path.join(config['dir']['log'],"subsample_k_{k}_PC_{pc}_round_{run_id}.log")
	params: jobname = "subsample_k_{k}_PC_{pc}_round_{run_id}",
		rate = config['stable']['subsample_rate'],
		resolution = config['stable']['subsample_resolution'].strip(),
		PreprocessSubsetData_pars = config.get("PreprocessSubsetData_subsample_pars", "")
	threads: 12
	message: "subsampling {params.rate} from the full data set, recluster using k of {wildcards.k} and {wildcards.pc} PCs for round {wildcards.run_id} using {threads} threads"
	script: "../scripts/clustStab/subsample.R"


## gather the subsampled and reclustered cell idents
rule gather_subsample:
	input: rds = SUBSAMPLE_K_PC
	output: path.join(config['dir']['data'],'evalClustStab',"gather_subsample.rds")
	log: path.join(config['dir']['log'],"gather_subsample.log")
	threads: 3
	message: "gathering idents for subsample k"
	script: "../scripts/clustStab/gather_subsample.R"

