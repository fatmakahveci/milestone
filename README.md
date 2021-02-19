### Table of Content

<!-- MarkdownTOC levels="1,2,3,4,5,6" autolink="true" bracket="round" autoanchor="true" style="ordered" indent="\t" -->

1. [Milestone Workflow](#milestone-workflow)
	1. [Installation](#installation)
	1. [Preprocess](#preprocess)
		1. [Downloading NCBI data to create the reference genome for given species](#downloading-ncbi-data-to-create-the-reference-genome-for-given-species)
		1. [Prodigal training file](#prodigal-training-file)
			1. [Downloading reference genome of given species to create prodigal training file](#downloading-reference-genome-of-given-species-to-create-prodigal-training-file)
	1. [Creation of MLST schema of the species](#creation-of-mlst-schema-of-the-species)
	1. [Reference genome creation](#reference-genome-creation)
	1. [Creation of MLST schema of the given sample](#creation-of-mlst-schema-of-the-given-sample)
		1. [1. SBG](#1-sbg)
		1. [2. VG](#2-vg)
		1. [Python Library Installation](#python-library-installation)

<!-- /MarkdownTOC -->

<a id="milestone-workflow"></a>
# Milestone Workflow

<div align="center"><img src="./images/milestone.png"></div>

- Milestone is an end-to-end workflow to create cgMLST schema of _given species_ using the set of assemblies provided by the user or NCBI's public database  and sample's raw reads.

<a id="installation"></a>
## Installation

- [ ] milestone
	+ `$ >> pip install milestone` 
- chewbbaca:
    + https://github.com/B-UMMI/chewBBACA
- snakemake:
    + `$ >> pip install snakemake`
- blastp (for chewBBACA):
    + `$ >> conda install -c bioconda blast`
    + add to `$PATH:/usr/local/ncbi/blast/bin`

<a id="preprocess"></a>
## Preprocess

<a id="downloading-ncbi-data-to-create-the-reference-genome-for-given-species"></a>
### Downloading NCBI data to create the reference genome for given species

- This step can be skipped if the user defined data set will be used.
- ncbi-genome-download: In case that public NCBI database is used for the creation of the MLST schema for the reference species.
    + `$ >> pip install ncbi-genome-download`
- Download genome assemblies of given species taxonomy ID:
    + `$ >> ncbi-genome-download -s refseq -l complete -T <tax_id> -F fasta bacteria`

<a id="prodigal-training-file"></a>
### Prodigal training file

- This step can be skipped if the training file is already available.
- Installation:
	+ [prodigal (github)](https://github.com/hyattpd/Prodigal/wiki/installation)
- Species reference genome download:
    + `$ >> bash download_species_reference_fasta.sh --species "<species_name>"`
- Create the training file via _prodigal_:
    + `$ >> prodigal -i <species_ncbi_name.fna> -t <species.trn> -p single`
- [ ] Available prodigal training files:
    + https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files

<a id="downloading-reference-genome-of-given-species-to-create-prodigal-training-file"></a>
#### Downloading reference genome of given species to create prodigal training file

- `$ >> bash download_species_reference_fasta.sh -s <species_name>`

<a id="creation-of-mlst-schema-of-the-species"></a>
## Creation of MLST schema of the species

- [ ] Available chewie MLST schemas:
    + https://chewbbaca.online/stats

- milestoneSnake.py

<a id="reference-genome-creation"></a>
## Reference genome creation

- [ ] Snakemake

- Convert gene and alleles FASTA files and sample set's cgMLST schema into reference VCF and FASTA files.

	```python
	 python3 create_reference_files.py --cglist ${samples_cgMLSTschema.txt} --fasta ${cg_allele_fasta_files} --cg_mlst_tsv ${samples_cgMLST.tsv} --strain_cluster ${strain_cluster_file.txt}
	```

<a id="creation-of-mlst-schema-of-the-given-sample"></a>
## Creation of MLST schema of the given sample

- [ ] Snakemake

<a id="1-sbg"></a>
### 1. SBG 

- alignment
	+ `$ >> time sbg-aligner -v similar_ref_2_samples.vcf --threads 8 -o ERR3464558.bam --reference similar_ref_2_samples.fasta -q ERR3464558_1.fastq -Q ERR3464558_2.fastq --read_group_library 'lib'` % 35 sec

- alignment quality check
	+ `$ >> samtools view -F 0x04 -f 0x2 -q 20 -b ERR3464558.bam | samtools sort -o ERR3464558.sorted.bam` % 1 min 10 sec

- remove PCR duplicates - > check freebayes.md (sam & bam input different)
	+ `$ >> samtools collate -o - ERR3464558.sorted.bam | samtools fixmate -m - - | samtools sort -o - - | samtools markdup -r - ERR3464558.sorted.rmdup.bam` % 1 min 34 sec

- index reference fasta and sample alignment bam
	+ `$ >> samtools faidx similar_ref_2_samples.fasta` % 1 sec
	+ `$ >> samtools index ERR3464558.sorted.rmdup.bam` % 3 sec

- bam -> fasta
	+ `$ >> freebayes -f similar_ref_2_samples.fasta ERR3464558.sorted.rmdup.bam | vcffilter -f "QUAL > 24" > ERR3464558.vcf` % 1 min 19 sec
	+ `$ >> bgzip ERR3464558.vcf && tabix -p vcf ERR3464558.vcf.gz` % 2 sec
	+ `$ >> cat similar_ref_2_samples.fasta | bcftools consensus ERR3464558.vcf.gz > ERR3464558.fasta` % 1 sec

<a id="2-vg"></a>
### 2. VG

- alignment
	+ `$ >> samtools faidx similar_ref_2_samples.fasta` % 1 sec
	+ `$ >> vg construct -r similar_ref_2_samples.fasta -v similar_ref_2_samples.vcf.gz > similar_ref_2_samples.vg`  % 1 sec
	+ `$ >> vg index -x similar_ref_2_samples.xg -g similar_ref_2_samples.gcsa -k 16 similar_ref_2_samples.vg` % 20 sec
	+ `$ >> vg map -x similar_ref_2_samples.xg -g similar_ref_2_samples.gcsa -f  ERR3464558_1.fastq -f ERR3464558_2.fastq --threads 8 --surject-to bam > ERR3464558.bam` % 1 hour 48 min 6 sec

- alignment quality check
	+ `$ >> samtools view -F 0x04 -f 0x2 -q 20 -b ERR3464558.bam | samtools sort -o ERR3464558.sorted.bam` % 1 min 9 sec

- remove PCR duplicates - > check freebayes.md (sam & bam input different)
	+ `$ >> samtools collate -o - ERR3464558.sorted.bam | samtools fixmate -m - - | samtools sort -o - - | samtools markdup -r - ERR3464558.sorted.rmdup.bam` % 1 min 32 sec

- index sample alignment bam
	+ `$ >> samtools index ERR3464558.sorted.rmdup.bam` % 2 sec

- bam -> vcf
	+ `$ >> freebayes -f similar_ref_2_samples.fasta ERR3464558.sorted.rmdup.bam | vcffilter -f "QUAL > 24" > ERR3464558.vcf` % 1 min 24 sec

- index vcf
	+ `$ >> bgzip ERR3464558.vcf && tabix -p vcf ERR3464558.vcf.gz` % 1 sec

- vcf -> fasta
	+ `$ >> cat similar_ref_2_samples.fasta | bcftools consensus ERR3464558.vcf.gz > ERR3464558.fasta`

<a id="python-library-installation"></a>
### Python Library Installation

> $ >> `pip3 install pipreqs`

> $ >> `pipreqs .`

> $ >> `pip3 install -r requirements.txt`

