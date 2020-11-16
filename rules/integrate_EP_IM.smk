if 'hard' in config['qc']['strategy']:
	inData='hardflted'
if 'outly' in config['qc']['strategy']:
	inData='outlyflted'
if 'pipe' in config['qc']['strategy']:
	inData='pipeCompflted'


rule merge_EP_IM:
	input:
		EP_list=path.join('config','samples_{cancer}_EP.txt'),
                IM_list=path.join('config','samples_{cancer}_IM.txt')
	output:
		path.join(config['dir']['data'],'{cancer}_'+inData+'.rds')
	params:
		indir=path.join(config['dir']['data'],'sctNormSeuObj_'+inData),
		infile="sctNormSeu_" + inData
	threads:
		20
	script:
		"../scripts/step7_.R"
		

rule integrate_all:
	input:
		expand(path.join(config['dir']['data'],'{cancer}_'+inData+'.rds'),cancer=['ESCC','GEJ'])
	output:
		path.join(config['dir']['log'],'merge_EP_IM.finish')
	shell:
		"touch {output}"


