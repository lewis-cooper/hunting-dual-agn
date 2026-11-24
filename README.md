# Close-separation quasar candidate search and vetting

Scripts for imaging candidate discovery, photometric QSO scoring, candidate vetting, and random-control tests. The imaging workflow starts from a Gaia-linked quasar master and local survey catalogues, finds the initial pairs, and generates its own component-association tables. The one-arcsecond route starts from the prepared 712-system input. Survey data are supplied separately; score/flag exports, cut tables, control caches and diagnostic figures are generated locally and excluded from Git.

The core selection uses either the magnitude-dependent or magnitude-independent model, propagates Quaia redshift uncertainties, and applies the paper's PM and DCR rules.

## Core scripts

| Task | Script |
| --- | --- |
| Find initial imaging pairs and build component associations | [`search_imaging_candidates.py`](scripts/search_imaging_candidates.py) |
| Score resolved companions | [`score_one_arcsecond_resolved.py`](scripts/score_one_arcsecond_resolved.py) |
| Score blended systems | [`score_one_arcsecond_unresolved.py`](scripts/score_one_arcsecond_unresolved.py) |
| Score imaging components | [`score_imaging_components.py`](scripts/score_imaging_components.py) |
| Calibrate the shared resolved-score threshold | [`calibrate_resolved_threshold.py`](scripts/calibrate_resolved_threshold.py) |
| Build one-arcsecond selections and cut tables | [`select_one_arcsecond_candidates.py`](scripts/select_one_arcsecond_candidates.py) |
| Build imaging selections and cut table | [`select_imaging_candidates.py`](scripts/select_imaging_candidates.py) |
| Crossmatch published systems and remove confirmed contaminants | [`crossmatch_literature.py`](scripts/crossmatch_literature.py) |
| Plot final-candidate separation versus redshift | [`plot_candidate_separation_redshift.py`](scripts/plot_candidate_separation_redshift.py) |

Other modules prepare local input tables, combine the selected searches and
produce score, colour and separation diagnostics. Shared model interfaces,
photometry handling and selection rules are in the root Python modules.

## Products

Running the workflow produces exactly six files in `products/`. There is one final candidate catalogue;
the two score audits contain every assessed input, including rejected systems.

| Product | Purpose |
| --- | --- |
| `final_candidates.csv` | The unique combined candidate catalogue after confirmed literature contaminant removal, including literature flags and citations. |
| `one_arcsecond_scores.csv` | Scores, route assignments and flags for all 712 one-arcsecond inputs. |
| `imaging_scores.csv` | Both-component scores and selection flags for the 1,600 imaging systems passed to scoring. |
| `one_arcsecond_cutflow.csv` | Cumulative counts for the resolved and unresolved routes. |
| `imaging_cutflow.csv` | Cumulative counts for the one-Gaia and two-Gaia imaging tiers. |
| `threshold_calibration.json` | The calibrated threshold and fit summary. |

These files are written under `PAPER1_OUTPUT_ROOT`, which defaults to
`products/`. Intermediate preparation tables, per-model calculations and plot
support data are written to a separate ignored working directory, defaulting
to `products_work/`. Set `PAPER1_WORK_ROOT` to override it. Working files are
not part of the distributed products and can be regenerated from the external
inputs. Route-specific selected-candidate lists are not written.

Read Gaia and candidate identifiers as strings. Missing scores have status or
reason fields and must not be interpreted as measured zero scores.

The reference input catalogues produce 443 unique systems before the literature stage:
100 resolved one-arcsecond, 44 unresolved one-arcsecond and 334 imaging systems,
with 35 shared systems. The bundled reference snapshot identifies 32 confirmed
contaminants, leaving **411 final candidates**: 88 resolved, 41 unresolved and
310 imaging systems, with 28 shared. Seven final candidates have verified
matches to confirmed dual AGN. Uncertain classifications are retained and flagged.

## Literature comparison and exclusions

Two curated reference CSVs are included under [`data/literature`](data/literature/README.md):
[`known_binaries.csv`](data/literature/known_binaries.csv) and
[`known_lenses_contaminants.csv`](data/literature/known_lenses_contaminants.csv).
They contain published positions, redshifts, separations, classifications and
publication provenance. Multiple image pairs or papers can describe the same
system; match counts deduplicate by candidate ID.

The automatic literature stage follows cross-pipeline deduplication. Confirmed
dual references use $0\leq z\leq4$ and $0<\theta\leq1.5''$. A positional finder
radius of $1''$ retrieves associations; redshift is not a matching criterion.
Pair identity requires two distinct Gaia IDs, agreement of both component
positions, or documented published pair geometry. A finder association alone
does not remove a candidate. The reference guide explains the astrometric
tolerances and qualified classifications, including the retained likely lens
J0823+2418.

Only verified, confirmed lenses, quasar–star/galaxy pairs, stellar pairs and
published projected quasar pairs are removed. Uncertain or conflicting
classifications remain flagged in the final catalogue. This is a literature
classification audit, not a measurement of sample completeness or purity.

To rerun this part after the selection stages:

```sh
python -m scripts.combine_candidate_catalogues
python -m scripts.crossmatch_literature
```

The single final catalogue stays in `products/final_candidates.csv`. Both
all-input score audits gain literature flags and `literature_final_pass`; their
original selection flags and scores retain their pre-literature meaning. Both
cut tables gain a final exclusion row, including both scoring-model variants.
The combined pre-literature table and detailed associations, candidate audit
and summary are kept in the ignored working directory under `PAPER1_WORK_ROOT`.
The literature stage runs offline from the two bundled CSVs.

## Selection rules

Resolved companions pass if an eligible model gives
$\ln\mathcal R\geq-6.81$. Imaging requires both components to pass under the
same model variant. Blended candidates require finite QQ/QS Bayes factors in
both variants and $\ln B_{\mathrm{QQ,QS}}>0$ in at least one.

The reference threshold is calibrated using 109 systems with finite score,
redshift and secondary $S_{\mathrm{PM}}<70$, before final PM selection.
The fitted mixture populations are calibration proxies rather than measured
catalogue purity or completeness.

All two-Gaia routes require secondary $S_{\mathrm{PM}}<3$ when finite. Missing
secondary PM is retained and flagged. Imaging also requires finite primary
$S_{\mathrm{PM}}<3$. The pair-DCR angle cut, finite and $>20^\circ$, is applied
once upstream to one-Gaia imaging systems. Two-Gaia systems are exempt.

The imaging selection retains $R_\perp\geq1.05$ using an **externally supplied
measurements table**. This release contains no image-profile measurement or
plotting code, cutout images or profile-building products. The required table
schema is documented in the [input guide](data/README.md).

[`redshift_uncertainty.py`](redshift_uncertainty.py) uses the published Quaia
per-object errors through `RedshiftMatch` for component scoring, retaining the
2000 km/s half-width and central-redshift support test. Blended Quaia systems
average QQ and QS likelihoods separately over nine Gaussian-weighted nodes,
clipped to [0.15, 4.35], retaining clipped duplicates and the finite-node
intersection. DESI and SDSS receive zero additional redshift uncertainty.

## Install dependencies

Use Python 3.11 or newer; the reference workflow was tested with Python 3.13.5.
From the repository root:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
git clone https://github.com/vasilybelokurov/qso_p_color.git ../qso_p_color
git -C ../qso_p_color checkout dbc392623d27e22744d5224f9f1a3fd52f6230a8
python -m pip install -e ../qso_p_color
export QSO_PCOLOR_PACKAGE_ROOT="$(cd ../qso_p_color && pwd)"
```

The pinned model package includes `20261005_magdep_qlf` and
`20261005_magindep_qlf` under `models/multisurvey_psf/`.

For scoring, download `SFD_dust_4096_ngp.fits` and `SFD_dust_4096_sgp.fits` from
the [SFD dust-map archive](https://doi.org/10.7910/DVN/EWCNL5), place both in an
`sfd/` subdirectory, and set `QSO_PCOLOR_DUSTMAPS_DIR` to its parent directory.

## Run with separately supplied data

The code can be inspected without the input catalogues. Executing the analysis
requires the local files listed in [data/README.md](data/README.md): the
prepared one-arcsecond inputs, the imaging master and survey measurements,
and the externally supplied morphology/visual-inspection decisions. The imaging
search is entirely offline. It reads local Tractor, Gaia and g-band CCD
catalogues with query-coverage metadata, then rebuilds the tables used by vetting.
The optional chance-alignment stage queries anonymous public Gaia/Legacy/MAST
services for its control sources.

```sh
export PAPER1_INPUT_ROOT=/path/to/prepared_inputs
export QSO_PCOLOR_DUSTMAPS_DIR=/path/to/dustmaps
export PAPER1_IMAGING_MORPHOLOGY_PATH=/path/to/imaging_morphology_measurements.csv
export PAPER1_VISUAL_REJECTIONS_PATH=/path/to/visual_inspection_rejections.csv
export PAPER1_OUTPUT_ROOT="$PWD/products_independent"
export PAPER1_FIGURE_ROOT="$PWD/figures_independent"
```

Copy [`imaging_search.example.json`](imaging_search.example.json) to
`imaging_search.local.json` and set its local catalogue/coverage paths. Paths
and glob patterns are resolved relative to the config file. Check availability
or run just the imaging search without loading the QSO models:

```sh
python -m scripts.search_imaging_candidates --config imaging_search.local.json --check-inputs
python -m scripts.run_analysis --stage imaging-search --resume
```

To find the initial pairs before Gaia/CCD data are available:

```sh
python -m scripts.search_imaging_candidates --config imaging_search.local.json --stage discover
```

This writes the discovered pairs and explicit `gaia_query_positions.csv` and
`ccd_query_positions.csv` request tables under `products_work/imaging_search/`.
It does not fetch those data. Complete Gaia/CCD inputs and their coverage logs
are required before the following association and vetting stages can run.

The search resolves overlapping DR9 north/south photometry, ranks the sources
within 1.5 arcsec of each Gaia parent, validates two positive-flux PSF components,
and matches each to local Gaia DR3 sources within 0.7 arcsec. The original fixed
primary assignment uses a redshift-local nearest-Tractor g-r reference with
the target itself excluded. DCR axes are rebuilt from quality-selected g-band
CCDs using their WCS, observing times and exposure/noise weights. These steps
preserve the existing search implementation; the later PM/DCR selection rules
remain in the vetting stages.

Generated imaging inputs live in `products_work/imaging_search/inputs/`, or
under the equivalent `PAPER1_WORK_ROOT`. Completed searches can be resumed only
when the source-file, code, parameter and generated-output hashes agree. The
master is supplied as an input catalogue; merging the original SDSS/DESI/Quaia
survey releases into that master is outside this workflow.

Individual stages run with `python -m scripts.MODULE` after their preceding
stages have produced the required intermediate files. To run the full local
search and vetting workflow after configuring its external inputs:

```sh
PAPER1_BLEND_WORKERS=4 sh run_frozen_analysis.sh
```

This runner checks external inputs, runs the imaging candidate search, validates
its generated tables, executes scoring/selection and literature exclusions,
and checks the paper's expected counts. It performs no imaging-data downloads.
Blended scoring memoizes identical population-prior evaluations in memory using
the same value-preserving cache as the control scorer. It still calculates every
system's likelihoods from the supplied photometry. Scoring is CPU intensive. Each
parallel worker loads its model assets. Set `PAPER1_PYTHON` to choose the
interpreter. For another sample, invoke the relevant modules and assess its
results independently of the paper-specific count checks.

The full workflow also plots the unique final candidates in projected proper
separation versus adopted redshift, with Planck18 curves for fixed angular
separations. To regenerate this figure directly from the final CSV:

```sh
python -m scripts.plot_candidate_separation_redshift
```

The plot is saved as `figures/joint/separation_redshift.pdf` and `.png` under
`PAPER1_FIGURE_ROOT`. It uses only `products/final_candidates.csv` under
`PAPER1_OUTPUT_ROOT`; shared systems are plotted once in their own group.
Use `--candidates /path/to/final_candidates.csv` and `--figure-dir /path/to/figures`
to override these paths. It needs no model assets or source catalogues.

Input catalogues, cutouts and generated outputs are excluded from Git.
The scripts generate the initial imaging pairs, score audits and final
catalogue locally from the configured inputs.

## Verification

Focused tests run without the paper's input catalogues:

```sh
python -m unittest discover -s verification -p 'test_*.py'
python verification/verify_redshift_uncertainty.py
```

With external input paths configured, check execution prerequisites using:

```sh
python verification/verify_input_contracts.py
```

After running the workflow with the paper's input catalogues and measurement
tables, its score outputs and figures can be checked with:

```sh
python verification/verify_example_run.py
```

The checks inspect the six products and figures. Detailed calibration-sample
and blended-likelihood checks also run when their working tables are available;
otherwise they are reported as skipped. No verification report is added to
`products/`. Testing supports implementation consistency; scientific
interpretation requires assessment of the method and source measurements.

## Run vetting followed by chance-alignment controls

The chance-alignment workflow uses the modules in
[`chance_alignments/engine`](chance_alignments/engine) for catalogue queries,
matching, scoring and statistical summaries. An input adapter connects these
calculations to the sequential vetting workflow.

The literature stage runs before the final-candidate plots. Chance-alignment
controls continue to use the input populations and scoring rules; they do not
use the literature-cleaned candidate sample to estimate method false positives.

After configuring the imaging search, vetting inputs and model assets above, copy
[`chance_alignment.example.json`](chance_alignment.example.json) to
`chance_alignment.local.json`. Set the paths to the full Gaia-linked quasar
master and published Quaia catalogue described in the [input guide](data/README.md).
Then run the complete sequence:

```sh
python -m scripts.run_analysis --chance-config chance_alignment.local.json --dry-run
python -m scripts.run_analysis --chance-config chance_alignment.local.json --resume
```

This runs the offline imaging search, vetting and literature stages first,
then constructs and scores chance-alignment controls. The latter can fetch
survey data; use `--stage vetting` for a local-only run. It prepares control
inputs from the generated vetting tables, samples nested 30,000/60,000 fields,
fetches survey controls, scores components and blends, and writes the statistical
report. A fresh run builds its control inputs from the vetting tables and
separately supplied parent catalogues; saved control data are optional.

To run just the control calculations after vetting:

```sh
python -m scripts.run_chance_alignments --config chance_alignment.local.json --resume
```

Individual control stages are `prepare`, `fetch`, `components`, `blends`, and
`report`, in that order. Use `--stage NAME` to run a stage separately. Logs,
queries and scoring checkpoints are kept under `chance_alignment_work/` by
default; the combined report is
`chance_alignment_work/run/paired_60000/RESULTS.md`. The six vetting products
remain in their existing output directory. Nothing in the control workspace
is committed to Git.

The controls measure **P(pass | unrelated matched control, known outcome)**.
They do not estimate contamination of a selected candidate catalogue. One-Gaia
imaging retains the empirical DCR survival factor from the real cut flow;
it is not a direct DCR-angle measurement on the re-paired controls. Missing
outcomes and repeated source identities are retained in the statistical analysis.
See the [control workflow guide](chance_alignments/README.md) for details.

For an offline replay using separately supplied saved control data:

```sh
python -m scripts.run_chance_alignments --replay-data /path/to/control_data --verify-reference
```

This recalculates the rates, uncertainty and figures without survey queries
or model evaluation. The saved control data are optional and excluded from Git.
