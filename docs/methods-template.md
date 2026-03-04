# Milestone Methods Template

Use this text as a starting point for papers, theses, or reports. Edit the bracketed placeholders and remove any claims that are not supported by your validation dataset.

## Example Methods Wording

Whole-genome multilocus sequence typing (wgMLST) profiles were generated with Milestone, a local schema-based allele-calling workflow for bacterial isolates. The schema consisted of [N] loci and was derived from [custom coding sequences / an imported public schema from PubMLST/BIGSdb / another source]. Schema provenance was recorded in `scheme_manifest.json`, and the schema was audited with Milestone's schema QC utilities before downstream analysis.

Raw paired-end reads were processed with Milestone's `allele_calling` workflow using the [aligner name] aligner and a minimum locus breadth threshold of [X]%. Locus calls reported as `LNF`, `ASM`, `ALM`, or other non-comparable states were treated as unresolved for downstream profile discrimination analyses. When a sample profile carried at least one differing comparable locus relative to another profile, the two isolates were considered distinguishable at the wgMLST profile level.

Pairwise profile comparisons and allele-distance summaries were produced with Milestone's `profile_compare` and `profile_matrix` utilities. These results were interpreted as wgMLST profile discrimination metrics rather than direct phylogenetic distance estimates.

## Reporting Checklist

- state whether the schema was custom, imported, or curated externally
- report the schema version or release identifier from `scheme_manifest.json`
- report the species and schema type (`wgmlst` or `cgmlst`)
- report the locus breadth threshold used for calling
- report how unresolved loci were handled
- avoid claiming that wgMLST discrimination alone proves evolutionary distance

## Claims To Avoid

- "Milestone identifies the exact strain regardless of relatedness"
- "wgMLST distance is equivalent to phylogenetic distance"
- "all novel allele calls represent biologically confirmed new alleles" without orthogonal validation
