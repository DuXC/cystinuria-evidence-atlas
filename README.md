# Cystinuria variant–structure–function evidence atlas

Data and analysis code supporting **Assay and allele context shape the functional evidence map for cystinuria variants**.

Xiancheng Du, Yongkun Zhu, Shuchun Tao, Ming Chen, Chunhui Liu and Chao Sun  
Department of Urology, Zhongda Hospital, School of Medicine, Southeast University, Nanjing, China

**Release:** [v1.0.0](https://github.com/DuXC/cystinuria-evidence-atlas/releases/tag/v1.0.0), 8 October 2026. The variant and structure snapshot is dated 5 October 2026; source curation extends through 8 October 2026.

## Contents and reuse

The atlas links 488 sequence-validated ClinVar VCV records with structural annotations and 105 selected result/condition rows from original studies. Six human protein experiment families supply 51 eligible protein result rows for 25 substitutions. Protein and transport evidence match 22 and 20 VCV records, respectively. Coverage is scoped to the selected sources.

| Location | Contents |
| --- | --- |
| [Data guide](supplementary_data/DATA_GUIDE.md) | Evidence units, joins, numeric fields and interpretation |
| [S1–S19](supplementary_data/) | Frozen supplementary tables with checksums and source lineage |
| [Analysis source data](analysis_source_data/) | Inputs for descriptive plots |
| [Methods](methods/SUPPLEMENTARY_METHODS.md) | Variant mapping, experiment eligibility and reuse examples |
| [Code](code/README.md) | Count verification and plotting instructions |
| [Figures](figures/README.md) | Six figures in PDF and editable SVG format |
| [References](references/references.tsv) | DOI-linked source bibliography |
| [Upstream sources](UPSTREAM_SOURCES.md) | Accessions, attribution and reuse conditions |

Start with a VCV record in S1, join its `gene` and `uniprot_change` to S2, then inspect `experimental_family_id`, `source_locator`, `numeric_kind` and assay conditions. RNA evidence additionally requires the exact nucleotide allele and transcript. Empty numeric cells indicate unreported or inapplicable values, not zero. Protein correspondence does not establish clinical pathogenicity; geometric contacts and AlphaMissense predictions remain separate annotation layers.

## Reproduce the descriptive results

```bash
python3 code/reproduce_release.py
```

This command uses Python's standard library and verifies all 19 table hashes. [Plotting instructions](code/README.md) specify the optional numerical and graphics dependencies. The package reproduces analyses from frozen derived tables; the source-interpretation and original genomic-mapping steps are documented in Methods.

## Cite the dataset

Du, X., Zhu, Y., Tao, S., Chen, M., Liu, C., Sun, C. (2026). Cystinuria variant–structure–function evidence atlas (v1.0.0) [Data set and code]. GitHub. https://github.com/DuXC/cystinuria-evidence-atlas/releases/tag/v1.0.0

Machine-readable author names, ORCIDs, date and version are supplied in [CITATION.cff](CITATION.cff). Cite the original experimental sources and reused databases when using their findings.

## Rights and attribution

Author-created curation, documentation and figures are shared under [CC BY 4.0](LICENSE-DATA.md); analysis code is shared under the [MIT License](code/LICENSE). Reused database annotations and quoted source wording retain their original rights and attribution as described in [UPSTREAM_SOURCES.md](UPSTREAM_SOURCES.md). DOI links and source locators identify the original articles. Source-file fields are provenance labels for the research archive.

## Funding

China Postdoctoral Science Foundation (2024M750457, Chunhui Liu); Jiangsu Provincial Research Project on Traditional Chinese Medicine and Integrated Chinese-Western Medicine (ZXFZ2026021, Chao Sun); Jiangsu Health International Exchange Program (2026, Chao Sun).
