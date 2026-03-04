# Milestone Limitations

Milestone is a research-oriented wgMLST workflow. The following constraints should be stated explicitly in academic and operational use.

## Scientific Limits

- Milestone performs schema-based allele calling and wgMLST profile discrimination; it is not a full phylogenetic inference system.
- A `different` profile comparison result means that two isolates are distinguishable by at least one comparable locus. It does not quantify evolutionary relatedness on its own.
- Thresholds such as minimum locus breadth must be validated per species, schema, and sequencing context.
- Novel allele calls should be interpreted cautiously unless they are confirmed by orthogonal analysis or external curation.

## Operational Limits

- Milestone does not currently provide the public nomenclature and hosted surveillance context of platforms such as EnteroBase.
- Milestone does not currently provide the full audit and governance surface of enterprise laboratory software.
- Imported public schemas are only as reliable as the source metadata and the local validation applied after import.

## Reporting Recommendation

When using Milestone in a manuscript or thesis, describe it as:

- a local wgMLST profiling workflow
- a schema-based allele-calling and profile-comparison system

Do not describe it as:

- a universal strain identification system
- a standalone phylogenetic analysis platform
