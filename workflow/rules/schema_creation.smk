###################################
## Author: @fatmakhv             ##
## The latest update: 17/10/2021 ##
## Aim: Schema creation          ##
###################################

ruleorder: create_reference_vcf_fasta > index_fasta > reheader_vcf > sort_filter_zip_index_vcf

rule create_reference_vcf_fasta:
    input:
        schema_dir = config["schema_dir"],
        code_dir = config["working_dir"]
    output:
        reference_vcf = config["reference_vcf"],
        reference_fasta = config["reference_fasta"]
    message: "Reference fasta, info.txt, and vcf files are being created..."
    threads: config["parameters"]["threads"]
    params:
        reference_info_txt = f'{config["reference"]}_info.txt',
        log_file = config["schema_creation_log_file"]
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo "Rule {rule} running create_reference.py on {input.schema_dir}." | tee -a {params.log_file}
        echo "Output files are {output.reference_vcf}, {params.reference_info_txt}, and {output.reference_fasta}. " | tee -a {params.log_file}
        python {input.code_dir}/scripts/create_reference.py --schema_dir {input.schema_dir} --reference_vcf {output.reference_vcf} --reference_fasta {output.reference_fasta} --reference_info {params.reference_info_txt} --threads {threads} 2>&1 | tee -a {params.log_file}
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        '''

rule index_fasta:
    input:
        reference_fasta = config["reference_fasta"]
    output:
        reference_fasta_fai = f'{config["reference_fasta"]}.fai'
    message:
        "Reference fasta is being indexed..."
    params:
        log_file = config["schema_creation_log_file"]
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo 'Rule {rule} indexing {input.reference_fasta} | tee -a {params.log_file}'
        samtools faidx {input.reference_fasta} --fai-idx {output.reference_fasta_fai}
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        '''

rule reheader_vcf:
    input:
        reference_fasta_fai = f'{config["reference_fasta"]}.fai'
    message:
        "Reference vcf file header is being fixed..."
    params:
        log_file = config["schema_creation_log_file"],
        reference_vcf = config["reference_vcf"],
        temp_reference_vcf = f'{config["reference_vcf"]}.'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo "Rule {rule} fixes the header of {params.reference_vcf} using {input.reference_fasta_fai}..." | tee -a {params.log_file}
        bcftools reheader -f {input.reference_fasta_fai} {params.reference_vcf} -o {params.temp_reference_vcf}
        mv {params.temp_reference_vcf} {params.reference_vcf}
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        '''

rule sort_filter_zip_index_vcf:
    input:
        reference_vcf = config["reference_vcf"]
    output:
        reference_vcf_gz = f'{config["reference_vcf"]}.gz',
        reference_vcf_gz_tbi = f'{config["reference_vcf"]}.gz.tbi'
    message: "VCF file is being compressed..."
    params:
        temp_reference_vcf = f'{config["reference_vcf"]}.',
        log_file = f'{config["schema_creation_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo "{input.reference_vcf} file is being compressed." | tee -a {params.log_file}
        echo "Output compressed file is {output.reference_vcf_gz}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        bcftools sort {input.reference_vcf} -Ov -o {params.temp_reference_vcf}
        mv {params.temp_reference_vcf} {input.reference_vcf}
        bgzip {input.reference_vcf} && tabix -p vcf {output.reference_vcf_gz}
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        '''
