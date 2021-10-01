<div align="left"> <h1> <img src="images/milestone.png" alt="milestone_logo"> MILESTONE </h1> </div>

Milestone is an end-to-end sample-based MLST profile creation workflow for bacterial species.

---

## Table of Contents

<!-- MarkdownTOC -->

- Tutorial
- Setting up the analysis
  - Creating conda environment \(optional\)
    - Installing conda environment \(optional\)
- Milestone Workflow
  - Schema Creation
  - Allele Calling
  - Parameters
    - 1. Snakemake \(common for both schema_creation and allele_calling\)
    - 2. Schema Creation
    - 3. Allele Calling
- Citation

<!-- /MarkdownTOC -->

---

## Tutorial

- Milestone works on raw reads. It doesn't require any assembly steps.

- This part will be extended.

---

## Setting up the analysis

- **Downloading data from NCBI (optional):** `bash download_species_reference_fasta.sh -s <species_name>`

### Creating conda environment (optional)

- **pip installation:** `sudo apt install python-pip`
- **dependencies:** _requirements.txt_, python 3.8.10
- **conda installation:** `pip install conda`

#### Installing conda environment (optional)

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create --name milestone bcftools=1.13 biopython=1.79 chewbbaca=2.7.0 htslib=1.13 fastp=0.12.4 freebayes=1.3.2 minimap2=2.17 pysam=0.16.0.1 samtools=1.13 snakemake=5.32.2 vg=1.34
```

- VG only have installation via conda for Linux.
- For macOS, you need to delete `vg=1.34` from the command above and install VG to your local by following the steps in https://github.com/vgteam/vg#building-on-macos as an additional step.

- You can activate the created environment to work in it:
  - `source activate milestone`

- When your analysis is done, you can deactivate the created environment:
  - `conda deactivate`
  - Your environment will be kept unless you remove it. You can use it again by activating with the line given above.

---

## Milestone Workflow

- Milestone has a fully-automated workflow.

![milestone workflow](images/milestone_ug_workflow.png)

### Schema Creation

`$ >> python milestone.py schema_creation [-h] [-n] [-p] [-s SNAKEFILE] [-t THREADS] [-F] [--ri] [--unlock] [-q] -r REFERENCE -o OUTPUT -g GENOME_DIR [-mt MLST_TYPE]`

- It creates its own `Snakefile`. There is no need to create it. After creating the required environment for milestone, you can directly use the commands without any installations.

- `-r REFERENCE` is the raw name of the reference file without any extensions, such as `.vcf`, `.fasta`, and `_info.txt`.


- It uses user-provided CDS sequences and their allele sequences, and creates reference files accordingly.
- `-o OUTPUT`
  - It is the directory to be created for the output `reference.vcf`, `reference.fasta`, and `reference_info.txt` files.
- Milestone creates reference-related files:

![Graph representation in files](images/graph_vcf.png)

- `reference.fasta` contains the reference sequences of each coding sequences (CDSs).
  - It is created at the beginning of the analysis using user-defined set of alleles and CDSs.
  - It doesn't require any update after sample analysis because it has already contains the reference CDS of the sample's identified allele.
- `reference.vcf` contains all variations found on each alleles. 
  - It is created at the beginning of the analysis using user-defined set of alleles and CDSs.
  - For each allele, the variations are detected by aligning the alleles to its reference CDS.
  - Optionally, `reference.vcf` can be updated in `allele_calling` for the further analysis to involve variations in sample's novel alleles.
- `reference_info.txt` contains all variations on each allele of each CDS. 
- ![allele to vcf](images/allele_to_vcf_github.png)
  - Position (POS), reference (REF), alternate (ALT), and quality (QUAL) information of each variation are separated by specific characters in each line, where each variation of each allele is separated by comma(`,`) given in the same line (`cdsName_alleleId`).
    - i.e. `cdsName_alleleId POS*REF>ALT-QUAL,POS*REF>ALT-QUAL`
      - Each comma-separated part `POS*REF>ALT-QUAL` represents a variation of an allele.
      - Each variation set on a line, `POS*REF>ALT-QUAL,POS*REF>ALT-QUAL` , represents an allele.
      - Each line represents a single allele of a single CDS.
  - Optionally, `reference_info.txt` can be updated via `--update_reference` parameter of `allele_calling` after each sample analysis to involve newly identified alleles of the analyzed sample for the further analysis.

### Allele Calling

`$ >> milestone.py allele_calling [-h] [-n] [-p] [-s SNAKEFILE] [-t THREADS] [-F] [--ri] [--unlock] [-q] -r REFERENCE -o OUTPUT [--aligner ALIGNER] -e READ1 -E READ2 [--ur]`

- Milestone assigns the allele ID for sample's sequence aligned to the CDS based on the following criteria:
  - **<ID_from_the_reference>** If there is a complete match between the variations of sample's aligned sequence to the CDS and the allele-defining variation set given in TEXT-formatted reference file, it assigns the allele ID equal to the matching allele ID in the reference file.
  - **LNF** If the depth of coverage of the sample's CDS is lower than the expected, it assigns LNF (Locus Not Found) as allele ID to the sample's allele.
  - **1** If the depth of coverage of the sample's aligned sequence is equal to and more than the expected value and the sample does not have any variations for the CDS locus, it assigns the allele ID equal to the reference's, which is the longest allele of the reference CDS.
  - If there is no match between the variations of sample's aligned sequence and the allele-defining variation set given in TEXT-formatted reference file, it checks the validity of the sample's aligned sequence to the CDS before declaring the sequence as a novel allele of the CDS.
    - **LNF** If the length of the sequence is not a multiplier of 3 and/or the aligned sequence to the CDS contains in-frame stop codon, invalid start codon, and invalid stop codon, it assigns allele ID as LNF because bacterial genomes do not contain exons and it is not a valid coding sequence.
    - **ASM** If the sequence passes the validation steps, but its length is smaller than 20\% of the length of locus allele length mode, it assigns ASM (Alleles Smaller than Mode) to the sample's allele.
    - **ALM** If the sequence passes the validation steps, but its length is larger than 20\% of the length of locus allele length more, it assigns ALM (Alleles Larger than Mode) to the sample's allele.
- It creates its own `Snakefile`. There is no need to create it. After creating the required environment for milestone, you can directly use the commands without any installations.
- It uses the same output directory with `schema_creation`. It takes the reference files from the directory given with `-o OUTPUT` parameter.
- `-r REFERENCE` is the raw name of the reference file without any extensions, such as `.vcf`, `.fasta`, and `_info.txt`.
- `-e READ1` is the complete name of the read with its directory.
- `-E READ2` is the complete name of the read with its directory.
- If `--ur` (`--update_reference`) parameter is used to update the reference files for the further analysis, both `reference_info.txt` and `reference.vcf` will be updated to represent the novel alleles belonging to the analyzed sample.
- Reference update is described below:

![reference update](images/update_reference.png)

### Parameters

#### 1. Snakemake (common for both schema_creation and allele_calling)

-  `-n, --dryrun, --dry-run`

  Do not execute anything, and display what would be done. If you have a very large workflow, use `--dry-run --quiet` to just print a summary of the DAG of jobs.

- `-p, --printshellcmds`

  Print out the shell commands that will be executed.

- `-s SNAKEFILE, --snakefile SNAKEFILE` 

  The workflow definition in form of a snakefile. Usually, you should not need to specify this. By default, Snakemake will search for 'Snakefile','snakefile', 'workflow/Snakefile', 'workflow/snakefile' beneath the current working directory, in this order. Only if you definitely want a different layout, you need to use this parameter.

- `-t THREADS, --threads THREADS, --set-threads THREADS`

  Overwrite thread usage of rules. This allows to fine-tune workflow parallelization. In particular, this is helpful to target certain cluster nodes by e.g. shifting a rule to use more, or less threads than defined in the workflow. Thereby, THREADS has to be a positive integer, and RULE has to be the name of the rule.

- `-F, --forceall`

  Force the execution of the selected (or the first) rule and all rules it is dependent on regardless of already created output.

- `--ri, --rerun-incomplete`

  Re-run all jobs the output of which is recognized as incomplete.

- `--unlock`

  Remove a lock on the working directory.

- `-q, --quiet`

  Do not output any progress or rule information.

#### 2. Schema Creation

- `-h, --help`

  Show the instructions.

- `-r REFERENCE, --reference REFERENCE`

  Name of reference file to be given without extension and directory. (Both VCF and FASTA file name of the reference.)

- `-o OUTPUT, --output OUTPUT`

  Directory to be created for the output files

- `-sn SCHEMA_NAME_WITH_ITS_DIRECTORY, --schema_name vSCHEMA_NAME_WITH_ITS_DIRECTORY`

  User-provided gene directory name to create schema


#### 3. Allele Calling

- `-h, --help`

  Show the instructions.

- `-r REFERENCE, --reference REFERENCE`

  Name of reference file to be given without extension and directory. (Both VCF and FASTA file name of the reference.)

- `-o OUTPUT, --output OUTPUT`

  Directory to be created for the output files

- `--aligner ALIGNER`

  Allele Calling and Reference Update - Graph Aligner option, sbg or vg. (default: vg)

- `-e READ1, --read1 READ1`

  Allele Calling and Reference Update - Sample first read including its directory

- `-E READ2, --read2 READ2`

  Allele Calling and Reference Update - Sample second read including its directory

-  `--ur, --update_reference`

  Allele Calling and Reference Update - Update <reference_info.txt> and <reference.vcf> after the alignment of the given sample.

---

## Citation

@todo
