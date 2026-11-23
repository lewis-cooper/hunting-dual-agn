# Close-separation quasar candidate vetting

Scripts for photometric QSO scoring, candidate vetting, and random-control tests in the one-arcsecond and imaging searches. Input catalogues are supplied separately; score/flag exports, cut tables, control caches and diagnostic figures are generated locally and are excluded from Git.

The core selection uses either the magnitude-dependent or magnitude-independent model, propagates Quaia redshift uncertainties, and applies the paper's PM and DCR rules.

## Core scripts

| Task | Script |
| --- | --- |
| Score resolved companions | [`score_one_arcsecond_resolved.py`](scripts/score_one_arcsecond_resolved.py) |
| Score blended systems | [`score_one_arcsecond_unresolved.py`](scripts/score_one_arcsecond_unresolved.py) |
| Score imaging components | [`score_imaging_components.py`](scripts/score_imaging_components.py) |
| Calibrate the shared resolved-score threshold | [`calibrate_resolved_threshold.py`](scripts/calibrate_resolved_threshold.py) |
| Build one-arcsecond selections and cut tables | [`select_one_arcsecond_candidates.py`](scripts/select_one_arcsecond_candidates.py) |
| Build imaging selections and cut table | [`select_imaging_candidates.py`](scripts/select_imaging_candidates.py) |

Other modules prepare local input tables, combine the selected searches and
produce score, colour and separation diagnostics. Shared model interfaces,
photometry handling and selection rules are in the root Python modules.

## Products

Running the workflow produces exactly six files in `products/`. There is one final candidate catalogue;
the two score audits contain every assessed input, including rejected systems.

| Product | Purpose |
| --- | --- |
| `final_candidates.csv` | The unique combined final candidate catalogue. |
| `one_arcsecond_scores.csv` | Scores, route assignments and flags for all 712 one-arcsecond inputs. |
| `imaging_scores.csv` | Both-component scores and selection flags for all 1,600 imaging inputs. |
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

For the paper's input sample, the workflow produced 100 resolved one-arcsecond, 44 unresolved
one-arcsecond and 334 imaging systems, with 35 shared systems and **443 unique
candidates**. The manuscript's 415-system sample additionally applies 28
literature-contaminant exclusions outside this vetting workflow.

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
requires the prepared local files listed in [data/README.md](data/README.md).
The vetting scripts read those files; acquisition of the vetting inputs is
performed separately. The optional chance-alignment stage queries anonymous
public Gaia/Legacy/MAST services for its control sources.

```sh
export PAPER1_INPUT_ROOT=/path/to/prepared_inputs
export QSO_PCOLOR_DUSTMAPS_DIR=/path/to/dustmaps
export PAPER1_IMAGING_MORPHOLOGY_PATH=/path/to/imaging_morphology_measurements.csv
export PAPER1_VISUAL_REJECTIONS_PATH=/path/to/visual_inspection_rejections.csv
export PAPER1_OUTPUT_ROOT="$PWD/products_independent"
export PAPER1_FIGURE_ROOT="$PWD/figures_independent"
```

Individual stages run with `python -m scripts.MODULE` after their preceding
stages have produced the required intermediate files. To run the paper's full
workflow with its prepared inputs:

```sh
PAPER1_BLEND_WORKERS=4 sh run_frozen_analysis.sh
```

This runner checks input availability, executes scoring and selection, and
checks the paper's expected counts. Blended scoring is CPU intensive. Each
parallel worker loads its model assets. Set `PAPER1_PYTHON` to choose the
interpreter. For another sample, invoke the relevant modules and assess its
results independently of the paper-specific count checks.

Input catalogues, cutouts and generated outputs are excluded from Git.
The scripts generate the score audits and final catalogue locally from the
prepared inputs.

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

After running the workflow with the paper's prepared inputs, its score
outputs and figures can be checked with:

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

The literature crossmatch is maintained separately and is not included in this repository.

After configuring the vetting inputs and model assets above, copy
[`chance_alignment.example.json`](chance_alignment.example.json) to
`chance_alignment.local.json`. Set the paths to the full Gaia-linked quasar
master and published Quaia catalogue described in the [input guide](data/README.md).
Then run the complete sequence:

```sh
python -m scripts.run_analysis --chance-config chance_alignment.local.json --dry-run
python -m scripts.run_analysis --chance-config chance_alignment.local.json --resume
```

This runs the existing vetting checks and stages in their original order,
then constructs and scores chance-alignment controls. It prepares control
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
