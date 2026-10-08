# Supplementary data guide

The frozen variant/structure snapshot is dated 5 October 2026. Source curation
includes the Brauers2006 full text on 8 October 2026. TABLE_INDEX.tsv identifies
each input and output by SHA-256. Tables are UTF-8 tab-delimited text.

Units and joins
- S1: one VCV record. `variation_id` is the stable numeric key; `vcv` retains its version.
- S2: one selected result or assay condition; `evidence_id` is the row key.
- S3: one eligible original human protein experiment family.
- S4: the original 54 VCV records; two now have protein-matched experiments.
- S5: one source-family/substitution cell, including empty cells. `evidence_ids` link to S2.
- S6: one displayed gene–protein substitution, deduplicated within experiment families.
- S7: one established minigene haplotype context, not an isolated current nucleotide allele.
- S8/S11: source discrepancy and nomenclature records for S7.
- S9: one selected comparison study, with its inspected source and applicability.
- S10: comparator-set counts with explicitly named units.
- S12/S13: AGT1 other-gene/wild-type context and reported discrepancies.

Interpretation of core fields
- Protein joins use `gene` plus `uniprot_change` in S1/S4, `protein_variant` in S2,
  and `variant` in S5/S6. A shared protein substitution can map to multiple VCVs.
- RNA mapping additionally requires the exact nucleotide allele and transcript.
- `clinical_group` and `clinical_classification_raw` retain all-condition germline
  ClinVar aggregates. `has_cystinuria_rcv` describes a condition link, not a new classification.
- `has_human_primary_fulltext_function` means eligible human protein evidence,
  including purification/stability. The transport and numeric-summary flags are narrower.
- `source_clinically_observed` describes ascertainment, including control observations.
  `construct_origin` separately records how an experiment was designed.
- `counts_as_primary_human_protein` selects 51 result rows / 25 substitutions;
  `counts_as_primary_human_transport` selects 43 rows / 22 substitutions.
- `eligible_for_variant_protein_mapping` selects experiments with an exact current
  protein match. Protein mapping supplies neither nucleotide nor clinical equivalence.
- `numeric_kind` distinguishes exact summaries, rounded values, censored values,
  group ranges and qualitative observations. Blank numeric fields are unreported
  or inapplicable. They are not zero. `upper_bound_inclusive` controls the bound.
- `biological_replicates_n` and `technical_replicates_n` preserve different units.
  `transfection_to_assay_hours` is distinct from `uptake_duration_min`.
- S5 codes: N numeric summary; C strict upper bound; R textual group range;
  Q qualitative transport; P purification/stability; L protein detection/localization;
  S allele-specific splicing. Codes do not encode effect magnitude or pathogenicity.
- Geometry uses deposited wild-type all-heavy-atom distances in angstroms, taking
  the minimum across equivalent subunit copies. Plotted C-alpha positions are display markers.
- `in_*` comparator flags denote protein-list membership, separately from experiments.
- `source_locator` / `curation_source_locator` identify the original table, panel or page.
  `source_file` preserves the local research-archive path; originals remain in that archive.
- `human_feedback_kind` preserves actual feedback provenance: two signed source-review
  rows and 40 author-confirmed rows. The remaining 63 rows have no recorded human feedback.

## Source dependence and assay context

S14 removes one original protein family at a time; lost records describe coverage breadth.
S15 tabulates the current 4/5/6 angstrom rule; S19 contains its VCV membership.
S16 retains 12 source-summary rows for six paired substitutions. SD, SEM and strict bounds are different fields; no missing dispersion is filled.
S17 describes assay applicability, without assigning a clinical evidence-strength score.
S18 supplies the numerators and denominators of the descriptive coverage plots.

Quick use: locate a VCV in S1, use gene + uniprot_change to locate protein results in S2, then inspect experimental_family_id, source_locator, numeric_kind and conditions. RNA transfer additionally requires the stated nucleotide constraint. An empty match is scoped to this source set.

Reproduce derived counts with python3 code/reproduce_release.py from the unpacked package. See code/README.md for plotting instructions.
