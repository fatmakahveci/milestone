<div align="left"> <h1> <img src="images/milestone.png" alt="milestone_logo"> MILESTONE </h1> </div>

---

<!-- MarkdownTOC -->

- __Table of Contents__
	- Preprocess
		- Downloading NCBI data to create the reference genome for given species
		- Prodigal training file
			- Downloading reference genome of given species to create prodigal training file
	- Creation of reference genome
	- Creation of CDS FASTA file of the given sample
	- Creation of MLST schema of the given sample
		- extra

<!-- /MarkdownTOC -->

---

Milestone is an end-to-end cgMLST profile creation workflow for given bacterial species. It only uses available genome assemblies of the species provided by the user or NCBI's public database, and raw reads of the given sample.

- [ ] how to add milestone to pip: `$ >> ` 
- chewbbaca:
  + https://github.com/B-UMMI/chewBBACA
- snakemake:
  + `$ >> pip install snakemake`
- blastp (for chewBBACA):
  + `$ >> conda install -c bioconda blast`
  + add to `$PATH:/usr/local/ncbi/blast/bin`

## Preprocess

### Downloading NCBI data to create the reference genome for given species

- This step can be skipped if the user defined data set will be used.

- ncbi-genome-download: In case that public NCBI database is used for the creation of the MLST schema for the reference species.
  + `$ >> pip install ncbi-genome-download`

- Download genome assemblies of given species taxonomy ID:
  + `$ >> ncbi-genome-download -s refseq -l complete -T <tax_id> -F fasta bacteria`

### Prodigal training file

- This step can be skipped if the training file is already available.

- Installation: [prodigal (github)](https://github.com/hyattpd/Prodigal/wiki/installation)
- Species reference genome download:
  + `download_species_reference_fasta.sh --species "<species_name>"`
- Create the training file via _prodigal_:
  + `prodigal -i <species_ncbi_name.fna> -t <species.trn> -p single`
- Available prodigal training files:
  + https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files

#### Downloading reference genome of given species to create prodigal training file

- `$ bash download_species_reference_fasta.sh -s <species_name>`

## Creation of reference genome

- Available chewie MLST schemas: https://chewbbaca.online/stats

> python3 milestone.py chewbbaca -d <input_data_directory> -t <number_of_threads> --report -a <prodigal_training_file.trn> -g <reference_genome_assemblies_directory> -p -r <to_be_created_reference_file_name_without_extension> --snakefile Snakefile -F

## Creation of CDS FASTA file of the given sample

> python3 milestone.py alignment -d <input_data_directory> -t <number_of_threads> --report -p -r <to_be_created_reference_file_name_without_extension> -e <sample_1.fastq> -E <sample_2.fastq> --aligner <vg/sbg> --snakefile Snakefile -F

<a id="creation-of-mlst-schema-of-the-given-sample"></a>

## Creation of MLST schema of the given sample

> python3 comp_milestone_chewie.py -cm <results_alleles.tsv> -cs <directory_name_of_cds_files_created_by_chewbbaca> -mr <sample.fasta> -sid '<sample_id>' -o comp_results -fq1 <sample_1.fastq> -fq2 <sample_2.fastq> --idt 0.95 --mf 0.20 --lc 10 --vf 0.30 --sam

---
### extra

allele_name  <length> <isItmultipleof3> <start_codon?> <stop_codon>
length of alleles must be 3n
check the first 3 bases are start and last 3 bases are stop

=======

- __Table of Contents__
	- Preprocess
		- Downloading NCBI data to create the reference genome for given species
		- Prodigal training file
			- Downloading reference genome of given species to create prodigal training file
	- Creation of reference genome
	- Creation of CDS FASTA file of the given sample
	- Creation of MLST schema of the given sample
		- extra

<!-- /MarkdownTOC -->

---

Milestone is an end-to-end cgMLST profile creation workflow for given bacterial species. It only uses available genome assemblies of the species provided by the user or NCBI's public database, and raw reads of the given sample.

- [ ] how to add milestone to pip: `$ >> ` 
- chewbbaca:
  + https://github.com/B-UMMI/chewBBACA
- snakemake:
  + `$ >> pip install snakemake`
- blastp (for chewBBACA):
  + `$ >> conda install -c bioconda blast`
  + add to `$PATH:/usr/local/ncbi/blast/bin`

## Preprocess

### Downloading NCBI data to create the reference genome for given species

- This step can be skipped if the user defined data set will be used.

- ncbi-genome-download: In case that public NCBI database is used for the creation of the MLST schema for the reference species.
  + `$ >> pip install ncbi-genome-download`

- Download genome assemblies of given species taxonomy ID:
  + `$ >> ncbi-genome-download -s refseq -l complete -T <tax_id> -F fasta bacteria`

### Prodigal training file

- This step can be skipped if the training file is already available.

- Installation: [prodigal (github)](https://github.com/hyattpd/Prodigal/wiki/installation)
- Species reference genome download:
  + `download_species_reference_fasta.sh --species "<species_name>"`
- Create the training file via _prodigal_:
  + `prodigal -i <species_ncbi_name.fna> -t <species.trn> -p single`
- Available prodigal training files:
  + https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files

#### Downloading reference genome of given species to create prodigal training file

- `$ bash download_species_reference_fasta.sh -s <species_name>`

## Creation of reference genome

- Available chewie MLST schemas: https://chewbbaca.online/stats

> python3 milestone.py chewbbaca -d <input_data_directory> -t <number_of_threads> --report -a <prodigal_training_file.trn> -g <reference_genome_assemblies_directory> -p -r <to_be_created_reference_file_name_without_extension> --snakefile Snakefile -F

## Creation of CDS FASTA file of the given sample

> python3 milestone.py alignment -d <input_data_directory> -t <number_of_threads> --report -p -r <to_be_created_reference_file_name_without_extension> -e <sample_1.fastq> -E <sample_2.fastq> --aligner <vg/sbg> --snakefile Snakefile -F

<a id="creation-of-mlst-schema-of-the-given-sample"></a>

## Creation of MLST schema of the given sample

> python3 comp_milestone_chewie.py -cm <results_alleles.tsv> -cs <directory_name_of_cds_files_created_by_chewbbaca> -mr <sample.fasta> -sid '<sample_id>' -o comp_results -fq1 <sample_1.fastq> -fq2 <sample_2.fastq> --idt 0.95 --mf 0.20 --lc 10 --vf 0.30 --sam

---
### extra

allele_name  <length> <isItmultipleof3> <start_codon?> <stop_codon>
length of alleles must be 3n
check the first 3 bases are start and last 3 bases are stop
