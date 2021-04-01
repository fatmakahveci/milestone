ruleorder: create_wgmlst_schema > call_allele > create_cgmlst_schema > create_reference_vcf_fasta

rule create_wgmlst_schema:
    input:
        genome_dir = config["genome_dir"]
    output:
        data_dir = directory(f'{config["data_dir"]}/schema_seed')
    message: "chewBBACA is creating whole genome MLST (wgMLST) schema."
    params: log_file = f"{config['logs']}/chewbbaca.log"
    threads: config["parameters"]["threads"]
    shell:
        '''
        echo "python milestone.py chewbbaca is running..." | tee {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "chewBBACA CreateSchema is running on {input.genome_dir} with {threads} threads." | tee -a {params.log_file}
        echo "Output files on {output.data_dir}" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        chewBBACA.py CreateSchema -i {input.genome_dir} -o {output.data_dir} --cpu {threads} 2>&1 | tee -a {params.log_file}
        '''

rule call_allele:
    input:
        genome_dir = config["genome_dir"],
        schema_seed_dir = f'{config["data_dir"]}/schema_seed',
    output:
        allele_call_dir = directory(f'{config["data_dir"]}/allele_call')
    message: "chewBBACA is calling alleles."
    params: log_file = f"{config['logs']}/chewbbaca.log"
    threads: config["parameters"]["threads"]
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "chewBBACA AlleleCall is running on {input.genome_dir} and {input.schema_seed_dir}/schema_seed with {threads} threads." | tee -a {params.log_file}
        echo "Output files on {output.allele_call_dir}" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        chewBBACA.py AlleleCall -i {input.genome_dir} -g {input.schema_seed_dir}/schema_seed -o {output.allele_call_dir} --cpu {threads} 2>&1 | tee -a {params.log_file}
        '''

rule create_cgmlst_schema:
    input:
        allele_call_dir = f'{config["data_dir"]}/allele_call'
    output:
        cgmlst_dir = directory(config["cgmlst_dir"])
    message: "chewBBACA is creating core genome MLST (cgMLST) schema."
    params: log_file = f"{config['logs']}/chewbbaca.log"
    threads: config["parameters"]["threads"]
    shell:
        '''
        repeatedLoci=$(ls data/allele_call/result*/RepeatedLoci.txt)
        resultsAllelesTsv=$(ls {input.allele_call_dir}/result*/results_alleles.tsv)
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "chewBBACA ExtractCgMLST is running on results_alleles.tsv and RepeatedLoci.txt with 0.95 threshold." | tee -a {params.log_file}
        echo "Output file is {output.cgmlst_dir}" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        chewBBACA.py ExtractCgMLST -i "$resultsAllelesTsv" --r "$repeatedLoci" --t 0.95 -o {output.cgmlst_dir} 2>&1 | tee -a {params.log_file}
        '''

rule create_reference_vcf_fasta:
    input:
        schema_seed_dir = f'{config["data_dir"]}/schema_seed',
        cgmlst_dir = config["cgmlst_dir"]
    output:
        reference_vcf = "data/"+config["reference"]+".vcf",
        reference_fasta = "data/"+config["reference"]+".fasta"
    message: "Reference fasta and vcf files are being created..."
    threads: config["parameters"]["threads"]
    params: log_file = f"{config['logs']}/chewbbaca.log"
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "scripts/create_reference.py is runnning on {input.cgmlst_dir} and {input.schema_seed_dir}/schema_seed." | tee -a {params.log_file}
        echo "Output files '{output.reference_vcf}' and '{output.reference_fasta}' are created. " | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        python3 scripts/create_reference.py --cgmlst_dir {input.cgmlst_dir} --schema_seed_dir {input.schema_seed_dir}/schema_seed --reference_vcf {output.reference_vcf} --reference_fasta {output.reference_fasta} --threads {threads} 2>&1 | tee -a {params.log_file}
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        '''
