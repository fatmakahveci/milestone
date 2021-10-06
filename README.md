<div align="left"> <h1> <img src="images/milestone.png" alt="milestone_logo"> MILESTONE </h1> </div>

Milestone is an end-to-end sample-based MLST profile creation workflow for bacterial species.

---

## Table of Contents

<!-- MarkdownTOC -->

- Milestone Workflow
  - Schema Creation
  - Allele Calling
- Citation

<!-- /MarkdownTOC -->

---

## Milestone Workflow

- Milestone has a fully-automated workflow.

![milestone workflow](images/milestone_ug_workflow.png)

### Schema Creation


- Milestone creates reference-related files:

![Graph representation in files](images/graph_vcf.png)

![allele to vcf](images/allele_to_vcf_github.png)

**Details of \<reference\>_info.txt** Position (POS), reference (REF), alternate (ALT), and quality (QUAL) information of each variation are separated by specific characters in each line, where each variation of each allele is separated by comma(`,`) given in the same line (`cdsName_alleleId`).

- i.e. `cdsName_alleleId POS*REF>ALT-QUAL,POS*REF>ALT-QUAL`
  - Each comma-separated part `POS*REF>ALT-QUAL` represents a variation of an allele.
  - Each variation set on a line, `POS*REF>ALT-QUAL,POS*REF>ALT-QUAL` , represents an allele.
  - Each line represents a single allele of a single CDS.

### Allele Calling

- Milestone assigns the allele ID for sample's sequence aligned to the CDS based on the following criteria:
  - **<ID_from_the_reference>** If there is a complete match between the variations of sample's aligned sequence to the CDS and the allele-defining variation set given in TEXT-formatted reference file, it assigns the allele ID equal to the matching allele ID in the reference file.
  - **LNF** If the depth of coverage of the sample's CDS is lower than the expected, it assigns LNF (Locus Not Found) as allele ID to the sample's allele.
  - **1** If the depth of coverage of the sample's aligned sequence is equal to and more than the expected value and the sample does not have any variations for the CDS locus, it assigns the allele ID equal to the reference's, which is the longest allele of the reference CDS.
  - If there is no match between the variations of sample's aligned sequence and the allele-defining variation set given in TEXT-formatted reference file, it checks the validity of the sample's aligned sequence to the CDS before declaring the sequence as a novel allele of the CDS.
    - **LNF** If the length of the sequence is not a multiplier of 3 and/or the aligned sequence to the CDS contains in-frame stop codon, invalid start codon, and invalid stop codon, it assigns allele ID as LNF because bacterial genomes do not contain exons and it is not a valid coding sequence.
    - **ASM** If the sequence passes the validation steps, but its length is smaller than 20\% of the length of locus allele length mode, it assigns ASM (Alleles Smaller than Mode) to the sample's allele.
    - **ALM** If the sequence passes the validation steps, but its length is larger than 20\% of the length of locus allele length more, it assigns ALM (Alleles Larger than Mode) to the sample's allele.

- Reference update is described below:

![reference update](images/update_reference.png)

## Citation

@todo
