###################################
## Author: @fatmakhv             ##
## The latest update: 17/10/2021 ##
## Aim: Calling alleles          ##
###################################


#####################################
## VG                              ##
#####################################

rule vg_index:
    input:
        reference_vcf_gz = f'{config["reference_vcf_gz"]}',
        reference_fasta = config["reference_fasta"]
    output:
        reference_dist = f'{config["aligner_reference"]}.dist',
        reference_giraffe_gbz = f'{config["aligner_reference"]}.giraffe.gbz',
        reference_min = f'{config["aligner_reference"]}.min'
    threads: config["parameters"]["threads"]
    message: "vg autoindex and vg giraffe are running..."
    threads: config["parameters"]["threads"]
    params:
        log_file = config["allele_calling_log_file"],
        reference = config["reference"],
        aligner_reference = config["aligner_reference"]
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo "vg autoindex is running on {input.reference_vcf_gz} and {input.reference_fasta} with {threads} threads." | tee -a {params.log_file}
        echo "Output files are {output.reference_dist}, {output.reference_giraffe_gbz}, and {output.reference_min}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        gunzip {input.reference_vcf_gz}; awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "LC_ALL=C sort -k1,1 -k2,2n"}}' {params.reference}.vcf > {params.reference}.sorted.vcf
        bcftools norm --rm-dup all {params.reference}.sorted.vcf > {params.reference}.vcf
        rm {params.reference}.sorted.vcf
        bgzip -f {params.reference}.vcf && tabix -p vcf {params.reference}.vcf.gz
        vg autoindex --ref-fasta {input.reference_fasta} --vcf {input.reference_vcf_gz} -t {threads} --workflow map --prefix {params.aligner_reference}
        vg autoindex --ref-fasta {input.reference_fasta} --vcf {input.reference_vcf_gz} -t {threads} --workflow giraffe --prefix {params.aligner_reference}
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        '''

rule vg_giraffe:
    input:
        reference_dist = f'{config["aligner_reference"]}.dist',
        reference_giraffe_gbz = f'{config["aligner_reference"]}.giraffe.gbz',
        reference_min = f'{config["aligner_reference"]}.min',
        read1 = config["samples"]["sample1"],
        read2 = config["samples"]["sample2"],
    output:
        sample_bam = f'{config["output_dir"]}/vg/{config["sample"]}.bam'
    threads: config["parameters"]["threads"]
    message: "VG giraffe is running..."
    params:
        reference =config["aligner_reference"],
        log_file = f'{config["allele_calling_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo "vg giraffe is running on {input.reference_dist}, {input.reference_giraffe_gbz}, {input.reference_min}, {input.read1}, and {input.read2} with {threads} threads." | tee -a {params.log_file}
        echo "Output file is {output.sample_bam}" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        fastp -i {input.read1} -I {input.read2} -o {input.read1}.filtered -O {input.read2}.filtered -n 0
        mv {input.read1}.filtered {input.read1}
        mv {input.read2}.filtered {input.read2}
        vg giraffe -Z {input.reference_giraffe_gbz} -m {input.reference_min} -d {input.reference_dist} -f {input.read1} -f {input.read2} -x {params.reference}.xg -o BAM > {output.sample_bam}
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        '''

#####################################
## SBG                             ##
#####################################

rule sbg_graf:
    input:
        reference_fasta = config["reference_fasta"],
        reference_vcf_gz = f'{config["reference"]}.vcf.gz',
        read1 = config["samples"]["sample1"],
        read2 = config["samples"]["sample2"],
    output:
        sample_bam = f'{config["output_dir"]}/sbg/{config["sample"]}.bam'
    params: log_file = f'{config["allele_calling_log_file"]}'
    threads: config["parameters"]["threads"]
    message: "SBG GRAF is running..."
    shell: 
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo "{input.reference_vcf_gz} is being uncompressed for SBG GRAF." | tee -a {params.log_file}
        echo "SBG GRAF is running on {input.reference_fasta}, {input.reference_vcf_gz}, {input.read1}, and {input.read2} with {threads} threads." | tee -a {params.log_file}
        echo "Output file is {output.sample_bam}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        fastp -i {input.read1} -I {input.read2} -o {input.read1}.filtered -O {input.read2}.filtered -n 0
        mv {input.read1}.filtered {input.read1}
        mv {input.read2}.filtered {input.read2}
        sbg-aligner -v {input.reference_vcf_gz} --threads {threads} --reference {input.reference_fasta} -q {input.read1} -Q {input.read2} --read_group_library 'lib' -o {output.sample_bam}
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        '''

#####################################
## FOR BOTH SBG AND VG             ##
#####################################

rule samtools_commands:
    input:
        sample_bam =  f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.bam'
    output:
        sample_bam_bai = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.bam.bai',
        sample_depth = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.depth'
    message: "Samtools commands are running..."
    params:
        sample = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}',
        log_file = f'{config["allele_calling_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo "samtools view, sort, collate, fixmate, and markdup are running on {input.sample_bam}." | tee -a {params.log_file}
        echo "Output files are {output.sample_bam_bai} and {output.sample_depth}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        samtools view -F 0x04 -f 0x2 -q 20 -b {input.sample_bam} | samtools sort -o - | samtools collate -o - - | samtools fixmate -m - - | samtools sort -o - | samtools markdup -r - - > {params.sample}.sorted.bam
        mv {params.sample}.sorted.bam {params.sample}.bam
        samtools index {params.sample}.bam {output.sample_bam_bai}
        samtools coverage {params.sample}.bam -o {output.sample_depth}
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        '''

rule bam_to_vcf:
    input:
        sample_bam =  f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.bam',
        sample_bam_bai = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.bam.bai',
        reference_fasta = config["reference_fasta"]
    output:
        sample_vcf = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.vcf'
    message: "Sample's variants are being called..."
    params: log_file = f'{config["allele_calling_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo "freebayes is running on {input.sample_bam}, {input.sample_bam_bai}, and {input.reference_fasta}." | tee -a {params.log_file}
        echo "Output file is {output.sample_vcf}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        freebayes -f {input.reference_fasta} {input.sample_bam} | vcffilter -f "QUAL > 24" > {output.sample_vcf}
        echo "REFERENCE" > REFERENCE
        bcftools reheader -s REFERENCE -o {output.sample_vcf}. {output.sample_vcf}
        mv {output.sample_vcf}. {output.sample_vcf}
        rm REFERENCE
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file} 
        '''

rule bam_to_sam:
    input:
        sample_bam =  f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.bam'
    output:
        sample_sam =  f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.sam'
    message: "Sample bam is being converted into sample sam"
    params: log_file = f'{config["allele_calling_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo "samtools view is running on {input.sample_bam}" | tee -a {params.log_file}
        echo "Output file is {output.sample_sam}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        samtools view -h {input.sample_bam} > {output.sample_sam}
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        '''

rule vcf_to_sample_allele_info:
    input:
        reference_fasta = config["reference_fasta"],
        reference_info_txt = config["reference_info_txt"],
        schema_dir = config["schema_dir"],
        sample_depth = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.depth',
        sample_vcf = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.vcf',
        sample_sam =  f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.sam',
        code_dir = config["working_dir"]
    output:
        sample_mlst = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}_mlst.tsv'
    message: "File for allele defining variants for each CDS is being created."
    threads: config["parameters"]["threads"]
    params:
        log_file = f'{config["allele_calling_log_file"]}',
        reference = config["reference"],
        update_reference = f'{config["update_reference"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo "create_sample_info_txt.py is running on {input.reference_fasta}, {input.reference_info_txt}, {params.reference}.vcf, {input.sample_depth}, {input.sample_sam}, {input.schema_dir}, and {input.sample_vcf}." | tee -a {params.log_file}
        echo "Reference info and vcf is being updated: {params.update_reference}." | tee -a {params.log_file}
        echo "Output file is {output.sample_mlst}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        python {input.code_dir}/scripts/create_sample_info_txt.py --sample_vcf {input.sample_vcf} --reference_info {input.reference_info_txt} --reference_vcf {params.reference}.vcf --sample_depth {input.sample_depth} --reference_fasta {input.reference_fasta} --update_reference {params.update_reference} --schema_dir {input.schema_dir} --sample_sam {input.sample_sam} --threads {threads}
        now=$(date +"%T")
        echo "End: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        '''

