<div align="left"> <h1> <img src="images/milestone.png" alt="milestone_logo"> Milestone Tutorial </h1> </div>

This tutorial aims to create multilocus sequence typing (MLST) from the user-defined coding sequences and raw reads. Begin the tutorial by creating the environment for milestone run by following the instructions below.

---

## Table of Contents

<!-- MarkdownTOC -->

- 1. Setup
	- 1.1. Setting up the data for the tutorial
	- 1.2. Setting up the environment for the tutorial
		- Linux
			- i. Install pip \(Pip Installs Packages\) using APT \(Advanced Packaging Tool\)
			- ii. Install conda
			- iii. Create the conda environment
		- macOS
			- i. Install homebrew \(The Missing Package Manager\)
			- ii. Install pip \(Pip Installs Packages\) using homebrew
			- iii. Install conda
			- iv. Create the conda environment
	- 2.1. `milestone.py schema_creation`
		- a. From genome assemblies of species
		- b. From coding sequences
			- b.1. Only coding sequences are available in the initial set.
			- b.2. Coding sequences and their alleles are available in the initial set.
		- 2.1.1. Input files
		- 2.1.2. Parameters
			- 2.1.2.a. Milestone parameters
			- 2.1.2.b. Snakemake parameters \(*optional\)
		- 2.1.3. Output files
	- 2.2. `milestone.py allele_calling`
		- 2.2.1. Input files
		- 2.2.2. Parameters
			- 2.2.2.a. Milestone parameters
				- 2.2.2.a.1. VG
				- 2.2.2.a.2. SBG
			- 2.2.2.b. Snakemake parameters \(*optional\)
		- 2.2.3. Output files
			- 2.2.3.1. VG
			- 2.2.3.1. SBG

<!-- /MarkdownTOC -->

---

## 1. Setup

### 1.1. Setting up the data for the tutorial

- Create a directory *milestone_tutorial* for these exercises.

- Copy files from ... into *milestone_tutorial* directory.

### 1.2. Setting up the environment for the tutorial

---

#### Linux

##### i. Install [pip](https://packaging.python.org/tutorials/installing-packages/) (Pip Installs Packages) using APT (Advanced Packaging Tool)

- `sudo apt-get update`
- `sudo apt-get install python3-pip`

##### ii. Install [conda](https://conda.io/)

- Follow the instructions in [conda](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html)'s website.

##### iii. Create the conda environment

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda create --name milestone bcftools=1.13 biopython=1.79 chewbbaca=2.7.0 htslib=1.13 fastp=0.12 freebayes=1.3.2 minimap2=2.22 pysam=0.16.0.1 samtools=1.13 snakemake=5.32.2 vg=1.34
```

- You can activate the created environment to work in it:
  - `source activate milestone`
- When your analysis is done, you can deactivate the created environment:
  - `conda deactivate`
  - Your environment will be kept unless you remove it. You can use it again by activating with the line given above.

---

#### macOS

##### i. Install [homebrew](https://brew.sh) (The Missing Package Manager)

`/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"`

##### ii. Install [pip](https://packaging.python.org/tutorials/installing-packages/) (Pip Installs Packages) using homebrew

- `brew install python3.8`

##### iii. Install [conda](https://conda.io/)

- Follow the instructions in [conda](https://conda.io/projects/conda/en/latest/user-guide/install/macos.html)'s website.

##### iv. Create the conda environment

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda create --name milestone bcftools=1.13 biopython=1.79 chewbbaca=2.7.0 htslib=1.13 fastp=0.12 freebayes=1.3.2 minimap2=2.22 pysam=0.16.0.1 samtools=1.13 snakemake=5.32.2 
```

- VG only have conda installation for Linux so you need to install VG to your local by following the steps in [VG](https://github.com/vgteam/vg#building-on-macos)'s website as an additional step.
- You can activate the created environment to work in it:
  - `source activate milestone`
- When your analysis is done, you can deactivate the created environment:
  - `conda deactivate`
  - Your environment will be kept unless you remove it. You can use it again by activating with the line given above.

---

- Milestone runs in two modes:

1. `python milestone.py schema_creation`
2. `python milestone.py allele_calling`

---

### 2.1. `milestone.py schema_creation`

- Milestone creates `Snakefile` file so it doesn't require to use `--snakefile SNAKEFILE` parameter. Only if you definitely want a different layout, you need to use this parameter.
- Milestone creates `config.yaml` files so you should not create this file.

---

#### a. From genome assemblies of species

- You can use [chewBBACA](https://github.com/B-UMMI/chewBBACA) to call alleles using public or user-provided genome assemblies belonging to the species.

- If you prefer using public genome assemblies of the species of the interest, you can download the public data by running `download_species_reference_fasta.sh`  script with the command below:

  `bash download_species_reference_fasta.sh -s <species_name>`

---

#### b. From coding sequences

##### b.1. Only coding sequences are available in the initial set.

- It appends all the coding sequences to create `<reference.fasta>` file.
- It creates an empty `<reference_info.txt>` file for further analysis.
- It creates a `<reference.vcf>` file containing only a default header for further analysis.

##### b.2. Coding sequences and their alleles are available in the initial set.

- It appends all the coding sequences to create `<reference.fasta>` file.
- It creates a `<reference_info.txt>` file to identify the allele set of user-provided coding sequences for further analysis.
- It creates a `<reference.vcf>` file containing a default header and variations between alleles and their reference coding sequences for further analysis.

---

#### 2.1.1. Input files

```pseudocode
schema_name
|- CDS1.fasta
|- CDS2.fasta
|- ...
|- CDSn.fasta
```

---

#### 2.1.2. Parameters

##### 2.1.2.a. Milestone parameters

```pseudocode
--reference <reference>
		Name of reference file to be given without extension and directory. It will be used to name <reference>.fasta, <reference>.vcf, and <reference>_info.txt files. [required: True]

--threads <threads>
		Number of threads to be used in the workflow. [default: 1]

--schema_name <schema_name>
		Name of directory containing user-defined coding sequences and their alleles. [required: True]

--output <output_directory>
		Name of directory for the output reference files. [required: True]
```

##### 2.1.2.b. Snakemake parameters (*optional)

```pseudocode
--help
		Show this help message and exit.

--dryrun
		Display the commands without running. [default: False]
		
--unlock
		Removes possible locks on Snakemake. [default: False]

--quiet
		Do not output any progress. [default: False]

--rerun-incomplete
		Rerun all incomplete jobs. [default: False]

--forceall
		Run all rules independent of being already created output. [default: False]

--printshellcmds
		Prints the shell command to be executed. [default: False]
```

---

#### 2.1.3. Output files

```pseudocode
output_directory
|- <reference>.fasta
|- <reference>.vcf.gz
|- <reference>_info.txt
```

---

### 2.2. `milestone.py allele_calling`

---

#### 2.2.1. Input files

- `schema_name` will only be used to add novel allele sequences identified using the analyzed sample.
- Allele info for the reference is already available in `<reference>_info.txt` file.

```pseudocode
output_directory
|- <reference>.fasta
|- <reference>.vcf.gz
|- <reference>_info.txt

<sample_read1>

<sample_read2>

schema_name
|- CDS1.fasta
|- CDS2.fasta
|- ...
|- CDSn.fasta
```

---

#### 2.2.2. Parameters

##### 2.2.2.a. Milestone parameters

###### 2.2.2.a.1. VG

```pseudocode
--aligner vg

--read1 <sample_read1>
		First read of sample including its directory. [required: True]

--read2 <sample_read2>
		Second read of sample including its directory. [required: True]

--reference <reference>
		Name of reference file to be given without extension and directory. It will be used to name <reference>.fasta, <reference>.vcf, and <reference>_info.txt files. [required: True]

--threads <threads>
		Number of threads to be used in the workflow. [default: 1]

--schema_name <schema_name>
		Name of directory containing user-defined coding sequences and their alleles. [required: True]

--output <output_directory>
		Name of directory for the output reference files. [required: True]

--update_reference
		Updates <reference>_info.txt and <reference>.vcf files by adding the information coming from the analyzed sample. [default: False]
```

###### 2.2.2.a.2. SBG

```pseudocode
--aligner sbg

--read1 <sample_read1>
		First read of sample including its directory. [required: True]

--read2 <sample_read2>
		Second read of sample including its directory. [required: True]

--reference <reference>
		Name of reference file to be given without extension and directory. It will be used to name <reference>.fasta, <reference>.vcf, and <reference>_info.txt files. [required: True]

--threads <threads>
		Number of threads to be used in the workflow. [default: 1]

--schema_name <schema_name>
		Name of directory containing user-defined coding sequences and their alleles. [required: True]

--output <output_directory>
		Name of directory for the output reference files. [required: True]

--update_reference
		Updates <reference>_info.txt and <reference>.vcf files by adding the information coming from the analyzed sample. [default: False]
```

##### 2.2.2.b. Snakemake parameters (*optional)

```pseudocode
--help
		Show this help message and exit.

--dryrun
		Display the commands without running. [default: False]
		
--unlock
		Removes possible locks on Snakemake. [default: False]

--quiet
		Do not output any progress. [default: False]

--rerun-incomplete
		Rerun all incomplete jobs. [default: False]

--forceall
		Run all rules independent of being already created output. [default: False]

--printshellcmds
		Prints the shell command to be executed. [default: False]
```

---

#### 2.2.3. Output files

- If `--update_reference` parameter is used, `<reference>_info.txt` and `<reference.vcf>` will be updated. (Note: `<reference>.fasta` is not required to be updated because the references of coding sequences will remain the same.)
- Files in `vg` or `sbg` are created from the scratch , but the remaining files are the updated versions of the input files.

##### 2.2.3.1. VG

```pseudocode
output_directory
|- <reference>_info.txt
|- <reference>.fasta
|- <reference>.vcf
|- vg
|--- <sample>.vcf
|--- <sample>_mlst.tsv
|--- <sample>.bam
|--- <sample>.depth

schema_name
|- CDS1.fasta
|- CDS2.fasta
|- ...
|- CDSn.fasta
```

##### 2.2.3.1. SBG

```pseudocode
output_directory
|- <reference>_info.txt
|- <reference>.fasta
|- <reference>.vcf
|- sbg
|--- <sample>.vcf
|--- <sample>_mlst.tsv
|--- <sample>.bam
|--- <sample>.depth

schema_name
|- CDS1.fasta
|- CDS2.fasta
|- ...
|- CDSn.fasta
```

---
