# Milestone Workflow

<div align="center"><img src="./milestone.png"></div>

#### Python Library Installation

- create requirements.txt file for the project

> $ >> `pip3 install pipreqs`

> $ >> `pipreqs .`

- install required python libraries
> $ >> `pip3 install -r requirements.txt`

:dna: **Step 0. Preprocess if you don't have complete genome set for schema creation:**

- Download complete genome data for the species of interest from NCBI database automatically, using *run_chewBBACA.py*.

- If you don't know taxonomy ID of your species of interest, you can search on [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=).
  ```python
  python3 run_chewbbaca.py --taxid=${your_species_tax_id}
  ```

- If you don't have training file for prodigal run, you can download NCBI reference FASTA file and run the following code that will create `data.trn` file:
	```python
	python3 run_chewbbaca.py --fasta=${your_species_reference_fasta_file}
	```

:dna: **Step 1. run_chewBBACA.py**

- Run chewBBACA to create both species and sample set's core genome multilocus sequence typing (cgMLST) schemas.

	```python
	python3 run_chewbbaca.py --cpu=${value} --genome_samples_dir=${complete_genome_directory} --samples_dir=${samples_FASTA_files_directory}
	```
	
- For more information:

	```python
	python3 run_chewbbaca.py --help
	```

:dna: **Step 2. create_reference_files.py**

- Strains and some samples related to this strains are required for training. If you don't have, you can use [PHYLOViZ](https://online.phyloviz.net/index) to cluster your samples and write the samples to a text file (i.e. strain_clusters.txt) using the following format:

---
STRAIN1:SAMPLE1,SAMPLE3

STRAIN2:SAMPLE2,SAMPLE5

STRAIN3:SAMPLE4

---

- Convert gene and alleles FASTA files and sample set's cgMLST schema into reference VCF and FASTA files.

	```python
	 python3 create_reference_files.py --cglist ${samples_cgMLSTschema.txt} --fasta ${cg_allele_fasta_files} --cg_mlst_tsv ${samples_cgMLST.tsv} --strain_cluster ${strain_cluster_file.txt}
	```

- For more information:

	```python
	python3 create_reference_files.py --help
	```

:dna: **Step 3. ...**
