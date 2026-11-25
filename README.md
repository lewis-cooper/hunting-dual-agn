# Hunting Dual AGN

Discover close quasar pairs from a supplied quasar catalogue, score their
photometry, apply the selection rules and remove confirmed literature
contaminants. The main entry point is `run_analysis.py`. WSDB is the
default catalogue source; equivalent local survey extracts support offline runs.

The workflow calculates Gaia neighbour counts, discovers the one-arcsecond and
imaging samples, builds component photometry and resolution classes, calibrates
the score threshold and measures image morphology. The automatic final catalogue
precedes visual inspection, which is performed outside this repository.

## Repository layout

| Directory | Contents |
| --- | --- |
| `dual_agn/` | Shared acquisition, discovery, measurement, scoring and selection code |
| `run_analysis.py` | Main entry point for either pipeline or both |
| `run_chance_alignments.py` | Run or replay the independent control experiment |
| `scripts/gaia_neighbour/` | Gaia companion discovery, resolved/blended scoring and selection |
| `scripts/imaging/` | Imaging pair discovery, component scoring, morphology and selection |
| `scripts/shared/` | Catalogue acquisition, input mapping, literature matching and final outputs |
| `examples/` | Main acquisition and chance-alignment configuration templates |
| `examples/local/` | Templates for locally supplied survey extracts |
| `data/` | Input documentation and cited literature reference catalogues |
| `chance_alignments/` | Independent method false-positive calculations |
| `verification/` | Automated tests and input/output checks |

Generated products, working data and figures are ignored by Git. Quasar and survey
catalogues, images and credentials are supplied separately.

## Install

Use Python 3.11 or newer. From the repository root:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt -r requirements-wsdb.txt
git clone https://github.com/vasilybelokurov/qso_p_color.git ../qso_p_color
git -C ../qso_p_color checkout dbc392623d27e22744d5224f9f1a3fd52f6230a8
python -m pip install -e ../qso_p_color
export QSO_PCOLOR_PACKAGE_ROOT="$(cd ../qso_p_color && pwd)"
```

The pinned package contains both model bundles under `models/multisurvey_psf/`.
Download the two SFD maps from the [SFD archive](https://doi.org/10.7910/DVN/EWCNL5),
place them in an `sfd/` directory, and set `QSO_PCOLOR_DUSTMAPS_DIR` to its parent.
The workflow checks models and dust maps before acquiring candidate data.
Local-file users can omit `requirements-wsdb.txt`.

## Configure the main run

Copy the template into the repository root:

```sh
cp examples/wsdb_inputs.example.json wsdb_inputs.local.json
```

Set `quasar_catalogue.catalogue` and its column map. A generic catalogue needs a
unique object ID, RA/Dec in degrees and redshift. Specify the redshift type;
photometric redshifts require uncertainties. Column mappings accommodate other
catalogues without edits to Python scripts. Relative paths are resolved from
the copied configuration's directory. The [input guide](data/README.md) covers
mixed catalogues and native survey identifiers.

The main run retrieves Gaia DR3, Tractor DR9, exact DESI DR1 target photometry,
PS1 DR1, NSC DR2 and SkyMapper DR4 through WSDB. Ground-based g-band exposure
metadata for DCR come from public NOIRLab Data Lab `ls_dr9.ccds_annotated`.
Missing morphology images are downloaded from the Legacy Surveys DR9 r-band
cutout service when configured.

Resolved and imaging scoring use Tractor grz. Blended scoring also uses available
PS1 Kron, NSC native and SkyMapper Petrosian measurements. Parent-quasar
identifiers and adopted redshifts are retained from the supplied catalogue.

WSDB credentials come from `WSDB_USER` and `WSDB_PASSWORD`, or an external
PostgreSQL service/password file. Defaults are host `wsdb.ast.cam.ac.uk` and
database `wsdb`. Override them with `WSDB_HOST`, `WSDB_PORT`, `WSDB_DATABASE` or
`WSDB_NAME`. Standard PostgreSQL service/environment settings are supported.
Keep credentials outside JSON and source files.

Check database table/column access without downloading catalogue rows:

```sh
python -m scripts.shared.acquire_inputs --config wsdb_inputs.local.json --check-schema
```

This connects to WSDB; it does not check the separate public exposure/image
services. Failed queries stop acquisition and are never recorded as empty cones.

## Run to final candidates

```sh
python run_analysis.py
```

The default is candidate vetting with both pipelines. Select a pipeline with:

```sh
python run_analysis.py --pipeline gaia_neighbour
python run_analysis.py --pipeline imaging
python run_analysis.py --pipeline both
```

`gaia_neighbour` includes both the resolved and unresolved one-arcsecond routes.
An imaging-only run also discovers and scores the Gaia-neighbour calibration
sample to determine the shared threshold. It skips Gaia blended scoring and
candidate selection, and its final catalogue contains imaging candidates only.
A Gaia-neighbour-only run skips imaging acquisition, morphology and selection.

This acquires raw data, discovers both samples, builds measurements, scores
systems, calibrates the threshold, applies PM/DCR/morphology cuts, deduplicates
physical systems, applies the literature audit and makes diagnostic figures.
Completion prints `Final candidate count` and the final CSV path.

Acquisition and analysis share discovery working folders. The following stages
reuse discovery only after checking input, code and output hashes. Gaia/CCD
association still runs after those data become available. The initial Tractor
search is not repeated unnecessarily.

Use `--dry-run` to print the sequence. Add `--resume` to reuse verified raw-query
and discovery checkpoints. Scoring and downstream selection run again. Fresh
acquisition computes its own neighbour counts and candidate membership; no fixed
historical list, score or morphology table is required.

Keep a test run separate with:

```sh
export PAPER1_OUTPUT_ROOT=/path/to/test/products
export PAPER1_WORK_ROOT=/path/to/test/work
export PAPER1_FIGURE_ROOT=/path/to/test/figures
python run_analysis.py --stage vetting --wsdb-config /path/to/test/wsdb_inputs.local.json
```

Discovery can also be run alone:

```sh
python run_analysis.py --stage one-arcsecond-search
python run_analysis.py --stage imaging-search
```

The one-arcsecond-only command retrieves Gaia data without supplementary
photometry or morphology images.

## Local survey extracts

Copy the three templates from `examples/local/` into the repository root as
`quasar_catalogue.local.json`, `imaging_search.local.json` and
`measurement_inputs.local.json`. Supply raw catalogues, complete query coverage
and image locations, then run:

```sh
python run_analysis.py --stage vetting --data-source local \
  --quasar-catalogue-config quasar_catalogue.local.json \
  --imaging-search-config imaging_search.local.json \
  --measurement-config measurement_inputs.local.json
```

This computes Gaia associations and neighbour counts. Set
`morphology_images.download_missing` to `false` for a fully offline candidate run.
Exports can be shared separately where redistribution is permitted and reused
through this route. See the [input guide](data/README.md).

## Products and selection

The default combined run writes six files under `products/`, or `PAPER1_OUTPUT_ROOT`:

| File | Contents |
| --- | --- |
| `final_candidates.csv` | Unique candidates after confirmed literature contaminant removal |
| `one_arcsecond_scores.csv` | Scores and eligibility, selection and literature flags for every discovered one-arcsecond input |
| `imaging_scores.csv` | Component scores and selection/literature flags for assessed imaging systems |
| `one_arcsecond_cutflow.csv` | Resolved and unresolved cumulative cuts |
| `imaging_cutflow.csv` | One-Gaia and two-Gaia cumulative cuts |
| `threshold_calibration.json` | Threshold and calibration summary |

A single-pipeline run writes four files: the final catalogue, its pipeline's
score audit and cutflow, and the threshold calibration. Gaia calibration
intermediates for imaging-only runs remain in the working folder. Use separate
output, work and figure folders for different pipeline choices. The runner stops
if the output folder contains products from an inactive pipeline.

Intermediates and raw queries stay under `PAPER1_WORK_ROOT`; figures stay under
`PAPER1_FIGURE_ROOT`. Defaults and scientific settings are in
`dual_agn/settings.py`. Preserve source identifiers as strings and missing-value
flags when reading outputs.

Selection uses either model. Imaging requires both components to pass under the
same variant. Threshold calibration precedes final PM selection and retains its
secondary PM significance limit of 70. Two-Gaia systems require secondary PM
significance below 3 when finite; missing secondary PM is retained and flagged.
Imaging also requires finite primary PM significance below 3. DCR angle greater
than 20 degrees applies once to one-Gaia imaging; two-Gaia systems are exempt.
Morphology requires `R_perp >= 1.05`. Photometric redshift uncertainties propagate
through both scorers; spectroscopic redshifts receive zero additional uncertainty.

The cited CSVs in [data/literature](data/literature/README.md) support the automatic
literature audit. Confirmed dual comparisons retain `0 <= z <= 4`, separation
at most 1.5 arcsec and a 1-arcsec positional finder. Removal requires verified
pair identity and a confirmed contaminant classification. Uncertain or conflicting
classifications remain flagged. Visual inspection is not applied.

## Chance-alignment controls

Copy `examples/chance_alignment.example.json` to `chance_alignment.local.json`
and set its full master and published Quaia paths. Run candidate selection
followed by the independent control experiment with:

```sh
python run_analysis.py --stage all --chance-config chance_alignment.local.json
```

To run controls after candidate selection:

```sh
python run_chance_alignments.py --config chance_alignment.local.json --resume
```

Controls estimate method false-positive rates from input populations, rather
than contamination within the final selected catalogue. They use public control
downloads and a separate ignored workspace. Details and offline replay commands
are in [chance_alignments/README.md](chance_alignments/README.md).
The control experiment compares both pipelines and requires `--pipeline both`.

## Verification

```sh
python -m unittest discover -s verification -p 'test_*.py'
python verification/verify_redshift_uncertainty.py
python verification/verify_input_contracts.py --models-only
```

Tests use local fixtures and mocked acquisition. They do not query surveys or run
the full population. After analysis, check output consistency with:

```sh
python verification/verify_example_run.py
```

Pass the same `--pipeline gaia_neighbour` or `--pipeline imaging` option to this
check after a single-pipeline run. The main workflow supplies it automatically.

A complete survey run is required to verify an actual candidate count. Different
masters, survey releases and coverage can change membership and calibration.
