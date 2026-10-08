# Upstream sources and attribution

The variant and structure snapshot was assembled on 5 October 2026. Literature curation extends through 8 October 2026. Original findings retain their source identity; published classifications and computational predictions are distinct from the authors' curation.

| Source | Use in this release | Access and rights |
| --- | --- | --- |
| NCBI ClinVar | Versioned VCV identifiers, HGVS, reported classifications and review status | Individual links are in S1. See [ClinVar data-use policy](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/); retain attribution to ClinVar and its submitters. |
| UniProt (P82251, Q07837) | Protein accessions, sequence positions and topology | [UniProt license](https://www.uniprot.org/help/license), CC BY 4.0. Sequence/topology annotations are subsetted and joined to the atlas. |
| wwPDB / RCSB PDB 6LI9 | Human transporter assembly and author-derived geometric contact distances | [Structure 6LI9](https://www.rcsb.org/structure/6LI9); [PDB data policy](https://www.rcsb.org/pages/policies), CC0 for archive data. Cite the original structural study listed in the bibliography. |
| AlphaMissense, Cheng et al. (2023) | Reused `am_pathogenicity` and `am_class` annotations in S1/S4 | [Original article](https://doi.org/10.1126/science.adg7492); [prediction source and license](https://github.com/google-deepmind/alphamissense#alphamissense-predictions-license), CC BY 4.0. Copyright 2023 DeepMind Technologies Limited. Scores are subsetted and joined; predictions are not clinical classifications. |
| Original experimental and comparison studies | Source-located factual result summaries and selected short quotations | DOI/PMID and table/panel/page locators are supplied in S2 and S7–S13. Cite the relevant originals. Their text and article files retain publisher or author rights. |

`source_file`, `source_files`, `primary_source_file` and related fields identify the archived input used in curation. They are provenance labels rather than download paths. The source bibliography provides public identifiers. Original article PDFs, reviewer-return workbooks and author contact records are held separately from this release.

The distributed scripts operate on derived TSV tables. Raw database queries, original genomic validation, structure-distance calculations and article interpretation are not rerun by these scripts. Consult Methods and the original sources for those upstream steps.
