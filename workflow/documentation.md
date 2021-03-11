### Table of Content

<!-- MarkdownTOC levels="1,2,3,4,5,6" autolink="true" bracket="round" autoanchor="true" style="ordered" indent="\t" -->

1. [README.md](#readmemd)
1. [Milestone](#milestone)
    1. [Preprocess](#preprocess)
        1. [Downloading NCBI data to create the reference genome for given species](#downloading-ncbi-data-to-create-the-reference-genome-for-given-species)
        1. [Prodigal training file](#prodigal-training-file)
            1. [Downloading reference genome of given species to create prodigal training file](#downloading-reference-genome-of-given-species-to-create-prodigal-training-file)
    1. [Creation of MLST schema of the species](#creation-of-mlst-schema-of-the-species)
    1. [Reference genome creation](#reference-genome-creation)
    1. [Creation of MLST schema of the given sample](#creation-of-mlst-schema-of-the-given-sample)
1. [Feb 08, 21](#feb-08-21)
    1. [FASTA OUTPUT](#fasta-output)
        1. [SBG](#sbg)
        1. [VG](#vg)
    1. [vcfstats](#vcfstats)
        1. [SBG](#sbg-1)
        1. [VG](#vg-1)
        1. [VG - with VG workflow](#vg---with-vg-workflow)
    1. [FASTA comparison w/ chewie](#fasta-comparison-w-chewie)
        1. [VG](#vg-2)
        1. [SBG](#sbg-2)

<!-- /MarkdownTOC -->

<a id="readmemd"></a>
# README.md

- [ ] how to add milestone to pip: `$ >> ` 
- chewbbaca:
    + https://github.com/B-UMMI/chewBBACA
- snakemake:
    + `$ >> pip install snakemake`
- blastp (for chewBBACA):
    + `$ >> conda install -c bioconda blast`
    + add to `$PATH:/usr/local/ncbi/blast/bin`

<a id="milestone"></a>
# Milestone

- Milestone is an end-to-end cgMLST profile creation workflow for given bacterial species and sample's raw reads.

- Milestone allows to create cgMLST schema of _given species_ for the set of assemblies provided by the user or NCBI's public database.

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

- Installation: [prodigal (github)](https://github.com/hyattpd/Prodigal/wiki/installation)
- Species reference genome download:
    + `download_species_reference_fasta.sh --species "<species_name>"`
- Create the training file via _prodigal_:
    + `prodigal -i <species_ncbi_name.fna> -t <species.trn> -p single`
- [ ] Ask Joao -> Available prodigal training files:
    + https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files

<a id="downloading-reference-genome-of-given-species-to-create-prodigal-training-file"></a>
#### Downloading reference genome of given species to create prodigal training file

- `$ bash download_species_reference_fasta.sh -s <species_name>`

<a id="creation-of-mlst-schema-of-the-species"></a>
## Creation of MLST schema of the species

- [ ] Ask Joao -> Available chewie MLST schemas:
    + https://chewbbaca.online/stats

<a id="reference-genome-creation"></a>
## Reference genome creation

- [ ] Snakemake

<a id="creation-of-mlst-schema-of-the-given-sample"></a>
## Creation of MLST schema of the given sample

- [ ] Snakemake

---

<a id="feb-08-21"></a>
# Feb 08, 21

<a id="fasta-output"></a>
## FASTA OUTPUT

- __lisa server__

<a id="sbg"></a>
### SBG 

1. chewie run (cgMLST schema and alleles)

2. species reference creation

3. alignment ======= DIFFERENT =======

- `$ >> time sbg-aligner -v similar_ref_2_samples.vcf --threads 8 -o ERR3464558.bam --reference similar_ref_2_samples.fasta -q ERR3464558_1.fastq -Q ERR3464558_2.fastq --read_group_library 'lib'` % 35 sec

4. alignment quality check

- `$ >> samtools view -F 0x04 -f 0x2 -q 20 -b ERR3464558.bam | samtools sort -o ERR3464558.sorted.bam` % 1 min 10 sec

5. remove PCR duplicates - > check freebayes.md (sam & bam input different)

- `$ >> samtools collate -o - ERR3464558.sorted.bam | samtools fixmate -m - - | samtools sort -o - - | samtools markdup -r - ERR3464558.sorted.rmdup.bam` % 1 min 34 sec

6. index reference fasta and sample alignment bam

- `$ >> samtools faidx similar_ref_2_samples.fasta` % 1 sec
- `$ >> samtools index ERR3464558.sorted.rmdup.bam` % 3 sec

7. bam -> fasta

- `$ >> freebayes -f similar_ref_2_samples.fasta ERR3464558.sorted.rmdup.bam | vcffilter -f "QUAL > 24" > ERR3464558.vcf` % 1 min 19 sec
- `$ >> bgzip ERR3464558.vcf && tabix -p vcf ERR3464558.vcf.gz` % 2 sec
- `$ >> cat similar_ref_2_samples.fasta | bcftools consensus ERR3464558.vcf.gz > ERR3464558.fasta` % 1 sec

---

<a id="vg"></a>
### VG

1. chewie run (cgMLST schema and alleles)

2. species reference creation

3. alignment ======= DIFFERENT =======

- `$ >> samtools faidx similar_ref_2_samples.fasta` % 1 sec

- `$ >> vg construct -r similar_ref_2_samples.fasta -v similar_ref_2_samples.vcf.gz > similar_ref_2_samples.vg`  % 1 sec

- `$ >> vg index -x similar_ref_2_samples.xg -g similar_ref_2_samples.gcsa -k 16 similar_ref_2_samples.vg` % 20 sec

- `$ >> vg map -x similar_ref_2_samples.xg -g similar_ref_2_samples.gcsa -f  ERR3464558_1.fastq -f ERR3464558_2.fastq --threads 8 --surject-to bam > ERR3464558.bam` % 1 hour 48 min 6 sec

4. alignment quality check

- `$ >> samtools view -F 0x04 -f 0x2 -q 20 -b ERR3464558.bam | samtools sort -o ERR3464558.sorted.bam` % 1 min 9 sec

5. remove PCR duplicates - > check freebayes.md (sam & bam input different)

- `$ >> samtools collate -o - ERR3464558.sorted.bam | samtools fixmate -m - - | samtools sort -o - - | samtools markdup -r - ERR3464558.sorted.rmdup.bam` % 1 min 32 sec

6. index sample alignment bam

- `$ >> samtools index ERR3464558.sorted.rmdup.bam` % 2 sec

7. bam -> vcf

- `$ >> freebayes -f similar_ref_2_samples.fasta ERR3464558.sorted.rmdup.bam | vcffilter -f "QUAL > 24" > ERR3464558.vcf` % 1 min 24 sec

8. index vcf

- `$ >> bgzip ERR3464558.vcf && tabix -p vcf ERR3464558.vcf.gz` % 1 sec

9. vcf -> fasta

- `$ >> cat similar_ref_2_samples.fasta | bcftools consensus ERR3464558.vcf.gz > ERR3464558.fasta`

<a id="vcfstats"></a>
## vcfstats

<a id="sbg-1"></a>
### SBG

total variant sites:    75
of which 74 (0.986667) are biallelic and 1 (0.0133333) are multiallelic
total variant alleles:  76
unique variant alleles: 76

snps:   66
mnps:   7
indels: 3
complex:    0

mismatches: 84

ts/tv ratio:    1.70968
deamination ratio:  0.827586
biallelic snps: 63 @ 1.73913

ins/del length frequency distribution
length  ins del ins/del
1   1   1   1
2           
3           
4           
5           
6           
7           
8           
9           
10          
11          
12          
13          
14          
15          
16          
17          
18          
19          
20          
21          
22          
23          
24          
25          
26          
27          
28          
29          
30          
31          
32          
33          
34          
35          
36      1   

insertion alleles / deletion alleles:   0.5
inserted bases / deleted bases: 0.027027

mnp length frequency distribution
length  count
2   1
3   3
4   
5   1
6   
7   2
total bases in mnps:    18

<a id="vg-1"></a>
### VG 

total variant sites:    92
of which 91 (0.98913) are biallelic and 1 (0.0108696) are multiallelic
total variant alleles:  93
unique variant alleles: 93

snps:   74
mnps:   17
indels: 2
complex:    0

mismatches: 119

ts/tv ratio:    1.76744
deamination ratio:  0.767442
biallelic snps: 73 @ 1.80769

ins/del length frequency distribution
length  ins del ins/del
1   1   1   1

insertion alleles / deletion alleles:   1
inserted bases / deleted bases: 1

mnp length frequency distribution
length  count
2   1
3   3
4   3
5   6
6   1
7   2
8   1
total bases in mnps:    45

<a id="vg---with-vg-workflow"></a>
### VG - with VG workflow

total variant sites:    56
of which 56 (1) are biallelic and 0 (0) are multiallelic
total variant alleles:  56
unique variant alleles: 56

snps:   54
mnps:   1
indels: 1
complex:    0

mismatches: 56

ts/tv ratio:    2.73333
deamination ratio:  0.863636
biallelic snps: 54 @ 2.6

ins/del length frequency distribution
length  ins del ins/del
1           
2       1   

insertion alleles / deletion alleles:   0
inserted bases / deleted bases: 0

mnp length frequency distribution
length  count
2   1
total bases in mnps:    2

runtime_vg_sbg.py

<a id="fasta-comparison-w-chewie"></a>
## FASTA comparison w/ chewie

<a id="vg-2"></a>
### VG

`$ >> time python3 comp_milestone_chewie.py -cm masked_masked_results_alleles.tsv -cs similar_schema -mr ERR3464558.fasta -sid 'ERR3464558' -o comp_results -fq1 ERR3464558_1.fastq -fq2 ERR3464558_2.fastq --idt 0.95 --mf 0.20 --lc 10 --vf 0.30 --sam` % 46 sec

<a id="sbg-2"></a>
### SBG

`$ >> python3 comp_milestone_chewie.py -cm masked_masked_results_alleles.tsv -cs similar_schema -mr ERR3464558.fasta -sid 'ERR3464558' -o comp_results -fq1 ERR3464558_1.fastq -fq2 ERR3464558_2.fastq --idt 0.95 --mf 0.20 --lc 10 --vf 0.30 --sam` % 49 sec

milestone --make-graph --vg --allele-vcf mystrains.vcf  --out myvggraph
milestone --call-alleles --vg --graph myvcfgraph --input input.fastq

allele_name  <length> <isItmultipleof3> <start_codon?> <stop_codon>
bunlar gen oldugu icin, intron da olmadigindan alel uzunluklari 3n olacakmis
bir de ilk 3 karakter start mi, son 3 karakter stop mu
