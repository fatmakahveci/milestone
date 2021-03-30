ruleorder: create_wgmlst_schema > call_allele > create_cgmlst_schema > create_reference_vcf_fasta

rule create_wgmlst_schema:
	input:
		genome_dir = config["genome_dir"]
	output:
		data_dir = directory(f'{config["data_dir"]}/schema_seed')
	message: "chewBBACA is creating whole genome MLST (wgMLST) schema."
	log: config["logs"]+"chewbbaca.log"
	threads: config["parameters"]["threads"]
	shell: "chewBBACA.py CreateSchema -i {input.genome_dir} -o {output.data_dir} --cpu {threads}"

rule call_allele:
	input:
		genome_dir = config["genome_dir"],
		schema_seed_dir = f'{config["data_dir"]}/schema_seed',
	output:
		allele_call_dir = directory(f'{config["data_dir"]}/allele_call')
	message: "chewBBACA is calling alleles."
	log: config["logs"]+"chewbbaca.log"
	threads: config["parameters"]["threads"]
	shell: "chewBBACA.py AlleleCall -i {input.genome_dir} -g {input.schema_seed_dir}/schema_seed -o {output.allele_call_dir} --cpu {threads}"

rule create_cgmlst_schema:
	input:
		allele_call_dir = f'{config["data_dir"]}/allele_call'
	output:
		cgmlst_dir = directory(config["cgmlst_dir"])
	message: "chewBBACA is creating core genome MLST (cgMLST) schema."
	log: "logs/chewbbaca.log"
	threads: config["parameters"]["threads"]
	shell:
		'repeatedLoci=$(ls data/allele_call/result*/RepeatedLoci.txt); resultsAllelesTsv=$(ls {input.allele_call_dir}/result*/results_alleles.tsv); chewBBACA.py ExtractCgMLST -i "$resultsAllelesTsv" --r "$repeatedLoci" --t 0.95 -o {output.cgmlst_dir}'

rule create_reference_vcf_fasta:
	input:
		schema_seed_dir = f'{config["data_dir"]}/schema_seed',
		cgmlst_dir = config["cgmlst_dir"]
	output:
		reference_vcf = "data/"+config["reference"]+".vcf",
		reference_fasta = "data/"+config["reference"]+".fasta"
	message: "Reference is being prepared..."
	threads: config["parameters"]["threads"]
	log: config["logs"]+"reference.log"
	shell: "python3 scripts/create_reference.py --cgmlst_dir {input.cgmlst_dir}  --schema_seed_dir {input.schema_seed_dir}/schema_seed --reference_vcf {output.reference_vcf} --reference_fasta {output.reference_fasta}"
