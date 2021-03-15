## VG ##

rule index_fasta:
	input:
		reference_fasta = "data/"+config["reference"]+".fasta"
	output:
		reference_fasta_fai = "data/"+config["reference"]+".fasta.fai"
	message: "Reference fasta file is being indexed..."
	log: config["logs"]+"reference.log"
	conda: config["envs"]+"defaults.yaml"
	shell: "samtools faidx {input.reference_fasta} -o {output.reference_fasta_fai} 2> {log}"

rule compress_vcf:
	input:
		reference_vcf = "data/"+config["reference"]+".vcf"
	output:
		reference_vcf_gz = "data/"+config["reference"]+".vcf.gz"
	message: "VCF file is being compressed..."
	log: config["logs"]+"reference.log"
	conda: config["envs"]+"defaults.yaml"
	shell: "bgzip {input.reference_vcf} && tabix -p vcf {output.reference_vcf_gz} 2> {log}"

rule vg_construct:
	input:
		reference_fasta = "data/"+config["reference"]+".fasta",
		reference_vcf_gz = "data/"+config["reference"]+".vcf.gz"
	output:
		reference_vg = "data/vg/"+config["reference"]+".vg"
	message: "VG construct is running..."
	log: config["logs"]+"reference.log"
	# conda: config["envs"]+"defaults.yaml"
	shell: "vg construct -r {input.reference_fasta} -v {input.reference_vcf_gz} > {output.reference_vg} 2> {log}"

rule vg_index:
	input:
		reference_vg = "data/vg/"+config["reference"]+".vg"
	output:
		reference_xg = "data/vg/"+config["reference"]+".xg",
		reference_gcsa = "data/vg/"+config["reference"]+".gcsa"
	threads: config["parameters"]["threads"]
	message: "VG index is running..."
	log: config["logs"]+"reference.log"
	# conda: config["envs"]+"defaults.yaml"
	shell: "vg index -x {output.reference_xg} -g {output.reference_gcsa} -k {threads} {input.reference_vg} 2> {log}"

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
	log: config["logs"]+"reference.log"
	# conda: config["envs"]+"defaults.yaml"
	shell: "vg map -x {input.reference_xg} -g {input.reference_gcsa} -f {input.read1} -f {input.read2} --threads {threads} --surject-to bam > {output.sample_bam} 2> {log}"

## SBG ##

rule sbg_graf:
	input:
		reference_fasta = "data/"+config["reference"]+".fasta",
		reference_vcf = "data/"+config["reference"]+".vcf",
		read1 = "data/"+config["samples"]["sample1"],
		read2 = "data/"+config["samples"]["sample2"]
	output:
		sample_bam = "data/sbg/"+"".join(config["samples"]["sample1"].split('_1')[0])+".bam"
	threads: config["parameters"]["threads"]
	log: "../logs/sbg.log"
	message: "SBG GRAF is running..." 
	log: config["logs"]+"reference.log"
	conda: config["envs"]+"sbg-graf.yaml"
	shell: "sbg-aligner -v {input.reference_vcf} --threads {threads} --reference {input.reference_fasta} -q {input.read1} -Q {input.read2} --read_group_library 'lib' -o {output.sample_bam}"

## SAMPLE.FASTA FOR BOTH SBG AND VG ##
rule alignment_quality_check:
	input:
		sample_bam = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".bam"
	output:
		sample_sorted_bam = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".sorted.bam"
	message: "Sample bam file is being sorted..."
	log: config["logs"]+"reference.log"
	conda: config["envs"]+"defaults.yaml"
	shell: "samtools view -F 0x04 -f 0x2 -q 20 -b {input.sample_bam} | samtools sort -o {output.sample_sorted_bam} 2> {log}"

# remove PCR duplicates - > check freebayes.md (sam & bam input different)
rule remove_pcr_duplicates:
	input:
		sample_sorted_bam = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".sorted.bam"
	output:
		sample_sorted_rmdup_bam = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".sorted.rmdup.bam"
	message: "PCR duplicates are being removed..."
	log: config["logs"]+"reference.log"
	conda: config["envs"]+"defaults.yaml"
	shell: "samtools collate -o - {input.sample_sorted_bam} | samtools fixmate -m - - | samtools sort -o - - | samtools markdup -r - {output.sample_sorted_rmdup_bam} 2> {log}"

rule index_sample_bam:
	input:
		sample_sorted_rmdup_bam = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".sorted.rmdup.bam"
	output:
		sample_sorted_rmdup_bam_bai = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".sorted.rmdup.bam.bai"
	message: "Sample bam file is being indexed..."
	log: config["logs"]+"reference.log"
	conda: config["envs"]+"defaults.yaml"
	shell: "samtools index {input.sample_sorted_rmdup_bam} {output.sample_sorted_rmdup_bam_bai} 2> {log}"

rule bam_to_vcf:
	input:
		sample_sorted_rmdup_bam = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".sorted.rmdup.bam",
		reference_fasta = "data/"+config["reference"]+".fasta"
	output:
		sample_vcf = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".vcf"
	message: "Sample's variant are being called..."
	log: config["logs"]+"reference.log"
	conda: config["envs"]+"defaults.yaml"
	shell: 'freebayes -f {input.reference_fasta} {input.sample_sorted_rmdup_bam} | vcffilter -f "QUAL > 24" > {output.sample_vcf} 2> {log}'

rule index_vcf:
	input:
		sample_vcf = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".vcf"
	output:
		sample_vcf_gz = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".vcf.gz"
	message: "Sample's VCF file is being compressed and indexed..."
	log: config["logs"]+"reference.log"
	conda: config["envs"]+"defaults.yaml"
	shell: "bgzip {input.sample_vcf} && tabix -p vcf {output.sample_vcf_gz} 2> {log}"

rule vcf_to_fasta:
	input:
		reference_fasta = "data/"+config["reference"]+".fasta",
		sample_vcf_gz = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".vcf.gz"
	output:
		sample_fasta = "data/"+config["aligner"]+"/"+"".join(config["samples"]["sample1"].split('_1')[0])+".fasta"
	message: "Sample's VCF file is being converted into sample's FASTA..."
	log: config["logs"]+"reference.log"
	conda: config["envs"]+"defaults.yaml"
	shell: "cat {input.reference_fasta} | bcftools consensus {input.sample_vcf_gz} -o {output.sample_fasta} 2> {log}"