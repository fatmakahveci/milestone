#####################################
## author: @fatmakhv               ##
## the latest update: June 14, 2021 ##
#####################################

#####################################
## VG                              ##
#####################################

rule index_fasta:
    input:
        reference_fasta = config["reference_fasta"]
    output:
        reference_fasta_fai = f'{config["reference_fasta"]}.fai'
    message: "Reference fasta file is being indexed..."
    params: log_file = f'{config["mlst_log_file"]}'
    shell:
        '''
        echo "'python milestone.py mlst' is running..." | tee {params.log_file}
        now=$(date +"%T")
        echo "Start: $now" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "{input.reference_fasta} file is being indexed." | tee -a {params.log_file}
        echo "Output index file is {output.reference_fasta_fai}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        samtools faidx {input.reference_fasta} -o {output.reference_fasta_fai}
        '''

rule compress_vcf:
    input:
        reference_vcf = config["reference_vcf"]
    output:
        reference_vcf_gz = f'{config["reference_vcf"]}.gz'
    message: "VCF file is being compressed..."
    params: log_file = f'{config["mlst_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "{input.reference_vcf} file is being compressed." | tee -a {params.log_file}
        echo "Output compressed file is {output.reference_vcf_gz}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        bgzip {input.reference_vcf} && tabix -p vcf {output.reference_vcf_gz}
        '''

rule decompress_vcf:
    input:
        reference_vcf_gz = f'{config["reference_vcf"]}.gz'
    output:
        reference_vcf = config["reference_vcf"]
    message: "VCF file is being decompressed..."
    params: log_file = f'{config["mlst_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "{input.reference_vcf_gz} file is being decompressed." | tee -a {params.log_file}
        echo "Output decompressed file is {output.reference_vcf}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        gunzip -f {input.reference_vcf_gz}
        '''

rule vg_construct:
    input:
        reference_fasta = config["reference_fasta"],
        reference_vcf_gz = f'{config["reference_vcf"]}.gz'
    output:
        reference_vg = f'{config["aligner_reference"]}.vg'
    threads: config["parameters"]["threads"]
    message: "VG construct is running..."
    params: log_file = f'{config["mlst_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "vg construct is running on {input.reference_fasta} and {input.reference_vcf_gz} files with {threads} threads." | tee -a {params.log_file}
        echo "Output file is {output.reference_vg}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        vg construct -r {input.reference_fasta} -v {input.reference_vcf_gz} -t {threads} > {output.reference_vg}
        '''

rule vg_index:
    input:
        reference_vg = f'{config["aligner_reference"]}.vg'
    output:
        reference_xg = f'{config["aligner_reference"]}.xg',
        reference_gcsa = f'{config["aligner_reference"]}.gcsa'
    threads: config["parameters"]["threads"]
    message: "VG index is running..."
    params: log_file = f'{config["mlst_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "vg index is running on {input.reference_vg} with {threads} threads." | tee -a {params.log_file}
        echo "Output files are {output.reference_xg} and {output.reference_gcsa}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        vg index -x {output.reference_xg} -g {output.reference_gcsa} -k {threads} -t {threads} {input.reference_vg}
        '''

rule vg_map:
    input:
        reference_xg = f'{config["aligner_reference"]}.xg',
        reference_gcsa = f'{config["aligner_reference"]}.gcsa',
        read1 = config["samples"]["sample1"],
        read2 = config["samples"]["sample2"]
    output:
        sample_bam = f'{config["output_dir"]}/vg/{config["sample"]}.bam'
    threads: config["parameters"]["threads"]
    message: "VG map is running..."
    params: log_file = f'{config["mlst_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "vg map is running on {input.reference_xg}, {input.reference_gcsa}, {input.read1}, and {input.read2} with {threads} threads." | tee -a {params.log_file}
        echo "Output file is {output.sample_bam}" | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        vg map -x {input.reference_xg} -g {input.reference_gcsa} -f {input.read1} -f {input.read2} -t {threads} --surject-to bam > {output.sample_bam}
        '''

#####################################
## SBG                             ##
#####################################

rule sbg_graf:
    input:
        reference_fasta = config["reference_fasta"],
        reference_vcf_gz = f'{config["reference"]}.vcf.gz',
        read1 = config["samples"]["sample1"],
        read2 = config["samples"]["sample2"]
    output:
        sample_bam = f'{config["output_dir"]}/sbg/{config["sample"]}.bam'
    params: log_file = f'{config["mlst_log_file"]}'
    threads: config["parameters"]["threads"]
    message: "SBG GRAF is running..."
    shell: 
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "{input.reference_vcf_gz} is being uncompressed for SBG GRAF." | tee -a {params.log_file}
        echo "SBG GRAF is running on {input.reference_fasta}, {input.reference_vcf_gz}, {input.read1}, and {input.read2} with {threads} threads." | tee -a {params.log_file}
        echo "Output file is {output.sample_bam}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        sbg-aligner -v {input.reference_vcf_gz} --threads {threads} --reference {input.reference_fasta} -q {input.read1} -Q {input.read2} --read_group_library 'lib' -o {output.sample_bam}
        '''

#####################################
## FOR BOTH SBG AND VG             ##
#####################################

rule alignment_quality_check:
    input:
        sample_bam =  f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.bam'
    output:
        sample_sorted_bam = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.sorted.bam'
    message: "Sample bam file is being sorted..."
    params: log_file = f'{config["mlst_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "samtools view and samtools sort are running on {input.sample_bam}." | tee -a {params.log_file}
        echo "Output file is {output.sample_sorted_bam}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        samtools view -F 0x04 -f 0x2 -q 20 -b {input.sample_bam} | samtools sort -o {output.sample_sorted_bam}
        '''

rule remove_pcr_duplicates:
    input:
        sample_sorted_bam = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.sorted.bam'
    output:
        sample_sorted_rmdup_bam = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.sorted.rmdup.bam'
    message: "PCR duplicates are being removed..."
    params: log_file = f'{config["mlst_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "PCR duplicates are being removed..." | tee -a {params.log_file}
        echo "samtools collate, samtools fixmate, samtools sort, and samtools markdup are running on {input.sample_sorted_bam}." | tee -a {params.log_file}
        echo "Output file is {output.sample_sorted_rmdup_bam}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        samtools collate -o - {input.sample_sorted_bam} | samtools fixmate -m - - | samtools sort -o - - | samtools markdup -r - {output.sample_sorted_rmdup_bam}
        '''

rule bam_to_sam:
    input:
        sample_sorted_bam = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.sorted.bam'
    output:
        sample_sam = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.sam'
    message: "Sample's bam file is being converted into sam format."
    params: log_file = f'{config["mlst_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "Sample's bam file is being converted into sam format." | tee -a {params.log_file}
        echo "samtools view is running on {input.sample_sorted_bam}." | tee -a {params.log_file}
        echo "Output file is {output.sample_sam}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        samtools view -h {input.sample_sorted_bam} -o {output.sample_sam}
        '''

rule index_sample_bam:
    input:
        sample_sorted_rmdup_bam = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.sorted.rmdup.bam'
    output:
        sample_sorted_rmdup_bam_bai = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.sorted.rmdup.bam.bai'
    message: "Sample bam file is being indexed..."
    params: log_file = f'{config["mlst_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "{input.sample_sorted_rmdup_bam} is being indexed." | tee -a {params.log_file}
        echo "Output file is {output.sample_sorted_rmdup_bam_bai}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        samtools index {input.sample_sorted_rmdup_bam} {output.sample_sorted_rmdup_bam_bai}
        '''

rule bam_to_vcf:
    input:
        sample_sorted_rmdup_bam = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.sorted.rmdup.bam',
        reference_fasta = config["reference_fasta"]
    output:
        sample_vcf = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.vcf'
    message: "Sample's variant are being called..."
    params: log_file = f'{config["mlst_log_file"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "freebayes is running on {input.sample_sorted_rmdup_bam} and {input.reference_fasta}." | tee -a {params.log_file}
        echo "Output file is {output.sample_vcf}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        freebayes -f {input.reference_fasta} {input.sample_sorted_rmdup_bam} | vcffilter -f "QUAL > 24" > {output.sample_vcf}
        '''

rule vcf_to_sample_allele_info:
    input:
        sample_vcf = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.vcf',
        sample_sam = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}.sam',
        reference_info_txt = config["reference_info_txt"],
        reference_fasta = config["reference_fasta"],
        reference_vcf = config["reference_vcf"],
        code_dir = config["working_dir"]
    output:
        sample_mlst = f'{config["output_dir"]}/{config["aligner"]}/{config["sample"]}_mlst.tsv'
    message: "File for allele defining variants for each CDS is being created."
    params:
        log_file = f'{config["mlst_log_file"]}',
        update_reference = f'{config["update_reference"]}'
    shell:
        '''
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "create_allele_info.py is running on {input.sample_vcf}, {input.sample_sam}, {input.reference_info_txt}, and {input.reference_fasta}." | tee -a {params.log_file}
        echo "Output file is {output.sample_mlst}." | tee -a {params.log_file}
        echo "---------------------------------------" | tee -a {params.log_file}
        echo "{params.update_reference}"
        python {input.code_dir}/scripts/create_allele_info.py --vcf {input.sample_vcf} --info_txt {input.reference_info_txt} --sam {input.sample_sam} --reference_fasta {input.reference_fasta} --reference_vcf {input.reference_vcf} --sample_mlst {output.sample_mlst} --update_reference {params.update_reference}
        '''