# Reproduction

Run `python3 code/reproduce_release.py` from the unpacked release. Python 3 standard library is sufficient. The script checks 19 frozen table hashes and reproduces protein/transport coverage, leave-one-family-out coverage, contact-threshold membership and paired numerical summaries. A successful run prints a JSON report with `"status": "PASS"`.

For Figures 1, 3 and S3, install `requirements-figures.txt` and run `python3 code/plot_figures.py`. The script writes regenerated figures into `figures/` and a plotting provenance record under `reproduction/`. The release also supplies the validated PDF and SVG versions of Figures 2, S1 and S2; this script does not regenerate those three illustrations.

Reproduction starts from the supplied derived variant and experiment tables. It does not rerun genomic validation, coordinate-distance extraction from raw structures or interpretation of the original articles. Source accessions, snapshot dates and reuse conditions are documented in `UPSTREAM_SOURCES.md`; joins and interpretation are described in `supplementary_data/DATA_GUIDE.md` and `methods/SUPPLEMENTARY_METHODS.md`.
