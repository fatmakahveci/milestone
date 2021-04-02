## VG ##

rule index_fasta:
	input:
		reference_fasta = "data/"+config["reference"]+".fasta"
	output:
		reference_fasta_fai = "data/"+config["reference"]+".fasta.fai"
	message: "Reference fasta file is being indexed..."
	params: log_file = f"{config['logs']}/mlst.log"
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
		reference_vcf = "data/"+config["reference"]+".vcf"
	output:
		reference_vcf_gz = "data/"+config["reference"]+".vcf.gz"
	message: "VCF file is being compressed..."
	params: log_file = f"{config['logs']}/mlst.log"
	shell:
		'''
		echo "---------------------------------------" | tee -a {params.log_file}
		echo "{input.reference_vcf} file is being compressed." | tee -a {params.log_file}
		echo "Output compressed file is {output.reference_vcf_gz}." | tee -a {params.log_file}
		echo "---------------------------------------" | tee -a {params.log_file}
		bgzip {input.reference_vcf} && tabix -p vcf {output.reference_vcf_gz}
		'''

rule vg_construct:
	input:
		reference_fasta = "data/"+config["reference"]+".fasta",
		reference_vcf_gz = "data/"+config["reference"]+".vcf.gz"
	output:
		reference_vg = "data/vg/"+config["reference"]+".vg"
	threads: config["parameters"]["threads"]
	message: "VG construct is running..."
	params: log_file = f"{config['logs']}/mlst.log"
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
		reference_vg = "data/vg/"+config["reference"]+".vg"
	output:
		reference_xg = "data/vg/"+config["reference"]+".xg",
		reference_gcsa = "data/vg/"+config["reference"]+".gcsa"
	threads: config["parameters"]["threads"]
	message: "VG index is running..."
	params: log_file = f"{config['logs']}/mlst.log"
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
		reference_xg = "data/vg/"+config["reference"]+".xg",
		reference_gcsa = "data/vg/"+config["reference"]+".gcsa",
		read1 = "data/"+config["samples"]["sample1"],
		read2 = "data/"+config["samples"]["sample2"],
	output:
		sample_bam = "data/vg/"+"".join(config["samples"]["sample1"].split('_1')[0])+".bam"
	threads: config["parameters"]["threads"]
	message: "VG map is running..."
	params: log_file = f"{config['logs']}/mlst.log"
	shell:
		'''
		echo "---------------------------------------" | tee -a {params.log_file}
		echo "vg map is running on {input.reference_xg}, {input.reference_gcsa}, {input.read1}, and {input.read2} with {threads} threads." | tee -a {params.log_file}
		echo "Output file is {output.sample_bam}" | tee -a {params.log_file}
		echo "---------------------------------------" | tee -a {params.log_file}
		vg map -x {input.reference_xg} -g {input.reference_gcsa} -f {input.read1} -f {input.read2} -t {threads} --surject-to bam > {output.sample_bam}
		'''

## SBG ##

rule sbg_graf:
	input:
		reference_fasta = "data/"+config["reference"]+".fasta",
		reference_vcf_gz = "data/"+config["reference"]+".vcf.gz",
		read1 = "data/"+config["samples"]["sample1"],
		read2 = "data/"+config["samples"]["sample2"]
	output:
		sample_bam = "data/sbg/"+"".join(config["samples"]["sample1"].split('_1')[0])+".bam"
	params:
		reference_vcf = f'data/{config["reference"]}.vcf'
	threads: config["parameters"]["threads"]
	params: log_file = f"{config['logs']}/mlst.log"
	message: "SBG GRAF is running..."
	shell: 
		'''
		echo "---------------------------------------" | tee -a {params.log_file}
		echo "{input.reference_vcf_gz} is being uncompressed for SBG GRAF." | tee -a {params.log_file}
		echo "SBG GRAF is running on {input.reference_fasta}, {params.reference_vcf}, {input.read1}, and {input.read2} with {threads} threads." | tee -a {params.log_file}
		echo "Output file is {output.sample_bam}." | tee -a {params.log_file}
		echo "---------------------------------------" | tee -a {params.log_file}
		gunzip {input.reference_vcf_gz}
		sbg-aligner -v {params.reference_vcf} --threads {threads} --reference {input.reference_fasta} -q {input.read1} -Q {input.read2} --read_group_library 'lib' -o {output.sample_bam}
		gzip {input.reference_vcf}
		'''

## SAMPLE.FASTA FOR BOTH SBG AND VG ##
rule alignment_quality_check:
	input:
		sample_bam = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".bam"
	output:
		sample_sorted_bam = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".sorted.bam"
	message: "Sample bam file is being sorted..."
	params: log_file = f"{config['logs']}/mlst.log"
	shell:
		'''
		echo "---------------------------------------" | tee -a {params.log_file}
		echo "samtools view and samtools sort are running on {input.sample_bam}." | tee -a {params.log_file}
		echo "Output file is {output.sample_sorted_bam}." | tee -a {params.log_file}
		echo "---------------------------------------" | tee -a {params.log_file}
		samtools view -F 0x04 -f 0x2 -q 20 -b {input.sample_bam} | samtools sort -o {output.sample_sorted_bam}
		'''

# remove PCR duplicates - > check freebayes.md (sam & bam input different)
rule remove_pcr_duplicates:
	input:
		sample_sorted_bam = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".sorted.bam"
	output:
		sample_sorted_rmdup_bam = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".sorted.rmdup.bam"
	message: "PCR duplicates are being removed..."
	params: log_file = f"{config['logs']}/mlst.log"
	shell:
		'''
		echo "---------------------------------------" | tee -a {params.log_file}
		echo "PCR duplicates are being removed..." | tee -a {params.log_file}
		echo "samtools collate, samtools fixmate, samtools sort, and samtools markdup are running on {input.sample_sorted_bam}." | tee -a {params.log_file}
		echo "Output file is {output.sample_sorted_rmdup_bam}." | tee -a {params.log_file}
		echo "---------------------------------------" | tee -a {params.log_file}
		samtools collate -o - {input.sample_sorted_bam} | samtools fixmate -m - - | samtools sort -o - - | samtools markdup -r - {output.sample_sorted_rmdup_bam}
		'''

rule index_sample_bam:
	input:
		sample_sorted_rmdup_bam = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".sorted.rmdup.bam"
	output:
		sample_sorted_rmdup_bam_bai = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".sorted.rmdup.bam.bai"
	message: "Sample bam file is being indexed..."
	params: log_file = f"{config['logs']}/mlst.log"
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
		sample_sorted_rmdup_bam = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".sorted.rmdup.bam",
		reference_fasta = "data/"+config["reference"]+".fasta"
	output:
		sample_vcf = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".vcf"
	message: "Sample's variant are being called..."
	params: log_file = f"{config['logs']}/mlst.log"
	shell:
		'''
		echo "---------------------------------------" | tee -a {params.log_file}
		echo "freebayes is running on {input.sample_sorted_rmdup_bam} and {input.reference_fasta}." | tee -a {params.log_file}
		echo "Output file is {output.sample_vcf}." | tee -a {params.log_file}
		echo "---------------------------------------" | tee -a {params.log_file}
		freebayes -f {input.reference_fasta} {input.sample_sorted_rmdup_bam} | vcffilter -f "QUAL > 24" > {output.sample_vcf}
		'''

rule index_vcf:
	input:
		sample_vcf = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".vcf"
	output:
		sample_vcf_gz = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".vcf.gz"
	message: "Sample's VCF file is being compressed and indexed..."
	params: log_file = f"{config['logs']}/mlst.log"
	shell:
		'''
		echo "---------------------------------------" | tee -a {params.log_file}
		echo "bgzip and tabix are running." | tee -a {params.log_file}
		echo "{input.sample_vcf} file is being compressed and indexed." | tee -a {params.log_file}
		echo "---------------------------------------" | tee -a {params.log_file}
		bgzip {input.sample_vcf} && tabix -p vcf {output.sample_vcf_gz}
		'''

rule vcf_to_fasta:
	input:
		reference_fasta = "data/"+config["reference"]+".fasta",
		sample_vcf_gz = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".vcf.gz"
	output:
		sample_fasta = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".fasta"
	message: "Sample's VCF file is being converted into sample's FASTA..."
	params: log_file = f"{config['logs']}/mlst.log"
	shell:
		'''
		echo "---------------------------------------" | tee -a {params.log_file}
		echo "bcftools consensus is running on {input.reference_fasta} and {input.sample_vcf_gz}." | tee -a {params.log_file}
		echo "Output file is {output.sample_fasta}" | tee -a {params.log_file}
		echo "---------------------------------------" | tee -a {params.log_file}
		bcftools consensus --fasta-ref {input.reference_fasta} {input.sample_vcf_gz} -o {output.sample_fasta}
		now=$(date +"%T")
		echo "End: $now" | tee -a {params.log_file}
		'''

# rule fasta_to_mlst:
# 	input:
# 		sample_fasta = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".fasta",
# 		allele_call_dir = config["allele_call_dir"],
# 		schema_seed_dir = f'{config["data_dir"]}/schema_seed',
# 		read1 = "data/"+config["samples"]["sample1"],
# 		read2 = "data/"+config["samples"]["sample2"]
# 		sid = (config["samples"]["sample1"]).split('_')[0]
# 	output:
# 		sample_mlst = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".tsv"
# 	message: "Sample's FASTA file is being converted into sample's MLST schema..."
# 	params: log_file = f"{config['logs']}/sample_mlst.log"
# 	shell:
# 		'''
# 		echo "---------------------------------------" | tee -a {params.log_file}
# 		echo "" | tee -a {params.log_file}
# 		echo "---------------------------------------" | tee -a {params.log_file}
# 		resultsAllelesTsv=$(ls {input.allele_call_dir}/result*/results_alleles.tsv)
# 		python3 scripts/convert_fasta_into_mlst.py -cm $resultsAllelesTsv -cs {input.schema_seed_dir} -mr {input.sample_fasta} \
# 		-sid '{params.sid}' -o {output.sample_mlst} -fq1 {input.read1} -fq2 {input.read2}
# 		'''
