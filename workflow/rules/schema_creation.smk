###################################
## Author: @fatmakhv             ##
## The latest update: 01/10/2021 ##
## Aim: Schema creation          ##
###################################

ruleorder: create_reference_vcf_fasta > compress_vcf

rule create_reference_vcf_fasta:
    input:
        schema_dir = config["schema_dir"],
        code_dir = config["working_dir"]
    output:
        reference_vcf = config["reference_vcf"],
        reference_fasta = config["reference_fasta"]
    message: "Reference fasta and vcf files are being created..."
    threads: config["parameters"]["threads"]
    params:
        reference_info_txt = f'{config["reference"]}_info.txt',
        log_file = config["schema_creation_log_file"]
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "create_reference.py is running on{input.schema_dir}." | tee -a {params.log_file}
        echo "Output files '{output.reference_vcf}', and '{output.reference_fasta}' are created. " | tee -a {params.log_file}
        pwd
        echo "---------------------------------------" | tee -a {params.log_file}
        python {input.code_dir}/scripts/create_reference.py --schema_dir {input.schema_dir} --reference_vcf {output.reference_vcf} --reference_fasta {output.reference_fasta} --reference_info {params.reference_info_txt} --threads {threads} 2>&1 | tee -a {params.log_file}
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        '''

rule compress_vcf:
    input:
        reference_vcf = config["reference_vcf"]
    output:
        reference_vcf_gz = f'{config["reference_vcf"]}.gz',
        reference_vcf_gz_tbi = f'{config["reference_vcf"]}.gz.tbi'
    message: "VCF file is being compressed..."
    params: log_file = f'{config["schema_creation_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "{input.reference_vcf} file is being compressed." | tee -a {params.log_file}
        echo "Output compressed file is {output.reference_vcf_gz}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        bgzip {input.reference_vcf} && tabix -p vcf {output.reference_vcf_gz}
        '''
