# Candidate catalogue

The publication catalogue contains 436 unique systems: 201 in the main tier and 235 in the additional tier.

- `final_candidates.csv`: complete curated catalogue.
- `main_tier.csv` and `additional_tier.csv`: disjoint parts of that same catalogue.
- `columns.csv`: field definitions and units.

The main tier requires two distinct Gaia detections. The additional tier contains one-Gaia imaging selections. The original identifiers are retained in `original_system_name` and `candidate_id`; coordinate names are display labels. Read Gaia source IDs as text, not floating-point numbers.

This release applies 1 explicit post-selection exclusion(s) to the automatic catalogue. The decisions and reasons are supplied to `export_catalogue.py` through its exclusions file. Automatic scores, selection rules and cutflows are unchanged. In this repository, `products/final_candidates.csv` is the automatic, pre-visual-inspection result; use the files in this folder for the publication sample.

`euclid_coverage.csv` records footprint membership, VIS pixel availability and the supplied visual classification for each publication candidate. `euclid_counts.csv` gives the resulting counts. Footprint membership does not guarantee usable imaging, and a resolved image pair does not confirm two AGN.
