
configfile: "./config.yaml"
samples = ['simu']
data_loc = '../align/'
bin_loc = '../bin/'
c_gap = 250

rule all:
	input:
		expand(data_loc + '{sample}.Aligned.sortedByCoord.out.addtag.bam', sample = samples),
		'polyApipe/simu_counts.tab.gz',
		'scAPA/scAPA.raw.Rds',
		'scAPAtrap/peak_assigned',
		'SCAPTURE/simu.KeepCell.assigned',
		'Sierra/sierra.Rds',
		'MAAPER/result.rds',
		'scDapars/dapars2_res.txt',
		'Sierra/Sierra.Rds',
		'benchmarks.RData'

rule generate_bed:
	output:
		bin_loc + 'filtered.fake_pa.bed'
	params:
		rscript = config['rscript'],
		fake_sc = config['generate_faked_pa'],
		bed = config['bed'],
		gap = c_gap
	shell:
		"""
		{params.rscript} {params.fake_sc} {params.bed} {output} {params.gap}
		"""

rule simu_fq:
	input:
		bin_loc + 'filtered.fake_pa.bed'
	output:
		fq = data_loc + 'simu.R2.fq.gz',
		gt = bin_loc + 'ground_truth.txt'
	params:
		script = config['simu_p'],
		fa = config['fasta']
	shell:
		"""
		zsh -c '
		. $HOME/.zshrc 
		conda activate apa_simu
		python {params.script} {params.fa} {input} {output.fq} {output.gt}		
		conda deactivate'

		"""

rule align:
	input:
		fq = data_loc + '{sample}.R2.fq.gz'
	output:
		data_loc + '{sample}.Aligned.sortedByCoord.out.bam',
		data_loc + '{sample}.Log.final.out'
	params:
		star = config['star'],
		index = config['index'],
		outprefix = data_loc + 'simu.'
	shell:
		"""
		{params.star} \
		--runThreadN 12 \
		--genomeDir {params.index} \
		--outSAMtype BAM SortedByCoordinate \
		--outBAMcompression 9 \
		--limitBAMsortRAM 60000000000 \
		--readFilesCommand zcat \
		--readFilesIn {input.fq} \
		--outFileNamePrefix {params.outprefix} && samtools index {output[0]}

		"""


rule addbc:
	input:
		bam = data_loc + '{sample}.Aligned.sortedByCoord.out.bam'
	output:
		bam = data_loc + '{sample}.Aligned.sortedByCoord.out.addtag.bam'
	params:
		pscript = config['p_script']['add_tag']
	shell:
		"""
		python {params.pscript} {input.bam} && samtools index {output.bam}
		"""

rule MAAPER:
	input:
		bam = expand(data_loc + '{sample}.Aligned.sortedByCoord.out.addtag.bam', sample = samples)
	output:
		'MAAPER/result.rds'
	params:
		pas_annotation = config['MAPPER']['pas_annotation'],
		gtf =  config['MAPPER']['gtf'],
		mapper =  config['MAPPER']['script'],
		prefix = 'MAAPER',
		rscript = config['rscript']
	shell:
		"""

		{params.rscript} {params.mapper} {params.pas_annotation} {params.gtf} {params.prefix} {input.bam}

		"""

rule polyApipe:
	input:
		bam = expand(data_loc + '{sample}.Aligned.sortedByCoord.out.addtag.bam', sample = samples)
	output:
		'polyApipe/simu_counts.tab.gz'
	params:
		prefix = 'polyApipe/simu',
	shell:
		"""
		zsh -c '
		. $HOME/.zshrc 
		conda activate polyApipe_env
		/home/zhouran/data/proj/2021-0918-apa_evaluation/04.evaluation_polyA/polyApipe/code/polyApipe-master/polyApipe.py -i \
		{input.bam} -o {params.prefix}
		conda deactivate'

		"""

rule scapa:
	input:
		bam = expand(data_loc + '{sample}.Aligned.sortedByCoord.out.addtag.bam', sample = samples)
	output:
		'scAPA/scAPA.raw.Rds'
	params:
		rscript = config['rscript'],
		script = config['scapa']['script'],
		prefix = 'scAPA'
	shell:
		"""
		{params.rscript} {params.script} {input.bam} {params.prefix}
		"""

rule scapatrap:
	input:
		bam = expand(data_loc + '{sample}.Aligned.sortedByCoord.out.addtag.bam', sample = samples)
	output:
		'scAPAtrap/peak_assigned'
	params:
		script = config['scapatrap']['script'],
		rscript = config['rscript'],
		prefix = 'scAPAtrap'
	shell:
		"""
		{params.rscript} {params.script} {input.bam} {params.prefix}
		"""

rule scapture:
	input:
		bam = expand(data_loc + '{sample}.Aligned.sortedByCoord.out.addtag.bam', sample = samples)
	output:
		'SCAPTURE/simu.KeepCell.assigned'
	params:
		script = config['scapture']['script'],
		prefix = 'SCAPTURE'
	shell:
		"""
		bash {params.script} {input.bam} {params.prefix}
		"""

rule Sierra:
	input:
		bam = expand(data_loc + '{sample}.Aligned.sortedByCoord.out.addtag.bam', sample = samples)
	output:
		'Sierra/sierra.Rds'
	params:
		script = config['sierra']['script'],
		rscript = config['rscript'],
		gtf =  config['MAPPER']['gtf'],
		prefix = 'Sierra'
	shell:
		"""
		{params.rscript} {params.script} {input.bam} {params.gtf} {params.prefix}
		"""

rule dapars_calculate_cov:
	input:
		bam = expand(data_loc + '{sample}.Aligned.sortedByCoord.out.addtag.bam', sample = samples)
	output:
		"scDapars/{sample}.wig"
	params:
		genomeCov = config["soft"]["genomeCov"]
	shell:
		"""
		{params.genomeCov} -bga -ibam {input} > \
		{output}
		"""
rule dapars_prepare_threeutr:
	input:
		bin_loc + 'filtered.fake_pa.bed'
	output:
		'scDapars/darpas2.bed'
	shell:
		"""
		awk 'BEGIN {{OFS="\t"}} {{print $1,$2,$3,$4,0,$6}}' \
		{input} \
		|bedtools sort -i /dev/stdin > {output}

		"""


rule dapars_config:
	input:
		log_file = expand(data_loc + '{sample}.Log.final.out', sample = samples),
		wig_file = expand('scDapars/{sample}.wig', sample = samples),
		three_utr_annotated = rules.dapars_prepare_threeutr.output
	output:
		seq_depth_file = 'scDapars/seq_depth',
		run_config = 'scDapars/run_config'
	params:
		num_threads = 12,
		cov_threshold = config['scdapars']['cov_threshold'],
		chromsome_id = config['scdapars']['chrom'],
		prefix = 'scDapars/output'

	run:
		out_lst = []

		for file in input.log_file:
			with open(file) as file_hd:
				for line in file_hd.readlines():
					line = line.strip()
					if line.startswith("Uniquely mapped reads number"):
						line = line.split("\t")
						out_lst.append("{}\t{}".format(file, line[-1]))
		with open(output.seq_depth_file, "w") as fileout_hd:
			fileout_hd.write("\n".join(out_lst))
		with open(output.run_config, "w") as config_hd:
			config_line = [
				'Annotated_3UTR={}'.format(input.three_utr_annotated),
				'Aligned_Wig_files={}'.format(','.join(input.wig_file)),
				'Output_directory={}'.format(params.prefix),
				'Output_result_file={}'.format('chromosome_'),
				'Coverage_threshold={}'.format(params.cov_threshold),
				'Num_Threads={}'.format(params.num_threads),
				'sequencing_depth_file={}'.format(output.seq_depth_file),
				'chromsome_id={}'.format(params.chromsome_id)
			]
			config_hd.write(
				'\n'.join(config_line)
			)			
rule run_dapars:
	input:
		seq_depth_file = 'scDapars/seq_depth',
		run_config = 'scDapars/run_config'
	output:
		'scDapars/welldone'
	params:
		p_script=config['scdapars']['p_script']
	shell:
		"""
		zsh -c '
		. $HOME/.zshrc # if not loaded automatically
		conda activate dapars
		{params.p_script} {input.run_config} {output}
		conda deactivate'
		"""

rule merge_res:
	input:
		wd = 'scDapars/welldone',
		bam = expand(data_loc + '{sample}.Aligned.sortedByCoord.out.addtag.bam', sample = samples)
	output:
		'scDapars/dapars2_res.txt'
	params:
		prefix = 'scDapars/output'
	shell:
		"""
		head -n 1 {params.prefix}/chromosome_10.txt > {output} && cat {params.prefix}/*.txt | grep -v "^Gene" >> {output}
		"""
rule make_bed:
	input:
		bin_loc + 'filtered.fake_pa.bed'
	output:
		bin_loc + 'collapse.bed'
	params:
		collapse_sc = config['collapse_sc']
	shell:
		"""
		python {params.collapse_sc} {input} {output}
		"""

rule prepare_res:
	input:
		bin_loc + 'collapse.bed',
		'polyApipe/simu_counts.tab.gz',
		'scAPA/scAPA.raw.Rds',
		'scAPAtrap/peak_assigned',
		'SCAPTURE/simu.KeepCell.assigned',
		'Sierra/sierra.Rds',
		'MAAPER/result.rds',
		'scDapars/dapars2_res.txt'
	output:
		'Sierra/Sierra.Rds'
	params:
		rscript=config['rscript'],
		prepare_sc = config['prepare_sc'],
		groud_truth = bin_loc + 'ground_truth.txt'
	shell:
		"""
		{params.rscript} {params.prepare_sc} {input[0]} {params.groud_truth}
		"""
rule compare:
	input:
		'Sierra/Sierra.Rds'
	output:
		'benchmarks.RData'
	params:
		rscript = config['rscript'],
		bm_sc = config['bm_sc']
	shell:
		"""
		{params.rscript} {params.bm_sc}
		"""
