# Hunting Dual AGN

Discover close quasar pairs, score their photometry, apply the candidate selection,
remove confirmed literature contaminants and write a final catalogue. The main
entry point is `scripts/run_analysis.py`. WSDB is the default catalogue source;
users without access can supply equivalent local survey extracts.

The starting input is a quasar catalogue. The workflow calculates Gaia neighbour
counts, discovers both search samples, builds component photometry and resolution
classes, calibrates the score threshold and measures image morphology. It does
not need a prepared 712-system list, candidate scores or morphology measurements.
Counts are computed for the supplied population. The automatic output precedes
visual inspection, which is performed outside this workflow.

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

The pinned package contains the magnitude-dependent and magnitude-independent
model bundles under `models/multisurvey_psf/`. Download the two SFD maps from the
[SFD archive](https://doi.org/10.7910/DVN/EWCNL5), put them in an `sfd/` directory,
and set `QSO_PCOLOR_DUSTMAPS_DIR` to its parent directory. Model/dust prerequisites
are checked before candidate data acquisition.

Local-file users can omit `requirements-wsdb.txt`. No input catalogue, image,
credential or generated analysis product is included in Git.

## Configure the WSDB run

Copy `wsdb_inputs.example.json` to `wsdb_inputs.local.json`. This is the one
configuration needed for candidate acquisition and analysis. Set its
`quasar_catalogue.catalogue` path and column map. A generic input needs a unique
object ID, RA/Dec in degrees and redshift. Specify the redshift kind explicitly;
photometric redshifts require their uncertainties. Other column names and larger
catalogues use the same scripts. The [input guide](data/README.md) covers mixed
catalogues, JSON fields and exact native survey identifiers.

The example retrieves raw records from the following sources:

| Data | Source |
| --- | --- |
| Gaia positions, magnitudes and PM measurements | WSDB Gaia DR3 |
| Tractor component photometry | WSDB Legacy Surveys DR9 |
| DESI target fallback photometry | WSDB DESI DR1 |
| PS1, NSC and SkyMapper supplementary photometry | WSDB PS1 DR1, NSC DR2 and SkyMapper DR4 |
| VHS and AllWISE inventory photometry | WSDB VHS DR5 and AllWISE |
| SDSS imaging photometry for the inventory | Public SDSS DR17 PhotoObjAll |
| Ground-based exposure metadata for DCR | Public NOIRLab Data Lab `ls_dr9.ccds_annotated` |
| Missing r-band morphology images | Legacy Surveys DR9 cutout service |

These preserve the catalogue releases and photometry types used by the existing
analysis: Tractor grz, PS1 Kron, SkyMapper Petrosian and native NSC magnitudes.
SDSS/VHS/AllWISE inventory measurements are not additional blended-scoring bands.
Current blended scoring adds PS1, NSC and SkyMapper to the available Tractor
system photometry. SDSS DR16Q parent identifiers/redshifts are distinct from the
SDSS DR17 imaging measurements.

WSDB connection details come from `WSDB_USER` and `WSDB_PASSWORD`, or an external
PostgreSQL service/password file. The defaults are host `wsdb.ast.cam.ac.uk` and
database `wsdb`. Override them with `WSDB_HOST`, `WSDB_PORT`, `WSDB_DATABASE` or
`WSDB_NAME`. Standard PostgreSQL service/environment settings are also supported.
Keep credentials outside JSON and source files. Local configuration, raw exports,
cutouts and analysis outputs are ignored by Git.

Check WSDB table/column access yourself without exporting catalogue rows:

```sh
python -m scripts.query_wsdb_inputs --config wsdb_inputs.local.json --check-schema
```

This check connects to WSDB. It does not validate availability of the separate
public image and exposure services. A successful empty query records an empty
cone; a failed query stops the acquisition and is never treated as a non-detection.

## Run from the quasar catalogue to final candidates

```sh
python -m scripts.run_analysis --stage vetting
```

The sequence acquires the configured raw data, discovers both samples, builds
measurements, scores systems, calibrates the threshold, applies PM/DCR/morphology
cuts, deduplicates physical systems, applies the literature audit and makes the
diagnostic figures. At successful completion it prints `Final candidate count`
and the final CSV path. `sh run_frozen_analysis.sh` invokes the same sequence.

Use a separate location for a fresh test:

```sh
export PAPER1_OUTPUT_ROOT=/path/to/test/products
export PAPER1_WORK_ROOT=/path/to/test/work
export PAPER1_FIGURE_ROOT=/path/to/test/figures
python -m scripts.run_analysis --stage vetting --wsdb-config /path/to/test/wsdb_inputs.local.json
```

Use `--dry-run` to print the sequence without running it. Add `--resume` after an
interruption to reuse verified raw query checkpoints. WSDB queries are batched;
query specifications, source releases, target positions, radii, code and output
hashes determine reuse. Raw source rows from overlapping batches are deduplicated
by the measurement builders. Neighbour counts and selection are recalculated;
there is no forced historical candidate count. Omit `--resume` for fresh queries.

For discovery alone:

```sh
python -m scripts.run_analysis --stage one-arcsecond-search
python -m scripts.run_analysis --stage imaging-search
```

The one-arcsecond-only command queries Gaia and computes discovery counts without
requesting Tractor cutouts or supplementary photometry.

## Local-file alternative

Supply native Gaia, Tractor and multi-survey catalogues, complete query coverage,
and local image paths using the three local examples. The master adapter computes
Gaia matches and neighbour counts from the raw sources:

```sh
python -m scripts.run_analysis --stage vetting --data-source local \
  --quasar-catalogue-config quasar_catalogue.local.json \
  --imaging-search-config imaging_search.local.json \
  --measurement-config measurement_inputs.local.json
```

This route has no WSDB connection. Set `morphology_images.download_missing` to
`false` for a fully offline candidate run. Generated WSDB exports can be shared
separately where redistribution is permitted and reused through this route.
The [input guide](data/README.md) explains path relocation and coverage formats.

## Products and scientific rules

There are exactly six files in `products/`, or `PAPER1_OUTPUT_ROOT`:

| File | Contents |
| --- | --- |
| `final_candidates.csv` | Unique combined candidates after confirmed literature contaminant removal |
| `one_arcsecond_scores.csv` | Scores, eligibility, selection and literature flags for every discovered one-arcsecond input |
| `imaging_scores.csv` | Component scores and selection/literature flags for assessed imaging systems |
| `one_arcsecond_cutflow.csv` | Resolved and unresolved cumulative cuts |
| `imaging_cutflow.csv` | One-Gaia and two-Gaia cumulative cuts |
| `threshold_calibration.json` | Threshold and calibration summary |

Raw data, query checkpoints, native associations and other intermediates stay
under `PAPER1_WORK_ROOT`; figures stay under `PAPER1_FIGURE_ROOT`. There is one
final candidate list. Read source identifiers as strings. Missing measurements
retain status/quality flags and must not be interpreted as measured zero values.

The selection uses either photometric model. Imaging requires both components
to pass under the same variant. The resolved threshold is calibrated before
final PM selection, retaining the calibration-only secondary PM significance
limit of 70. Two-Gaia systems require secondary PM significance below 3 when
finite; missing secondary PM is retained and flagged. Imaging also requires
finite primary PM significance below 3. DCR angle greater than 20 degrees is
applied once to one-Gaia imaging systems; two-Gaia systems are exempt.

Morphology uses the existing r-band profile-moment calculation and requires
`R_perp >= 1.05`. Missing images are downloaded only when configured. Their
coordinates/WCS are checked before measurement, and image provenance is retained.
Quaia errors propagate through component and blended scoring; spectroscopic
redshifts receive zero additional redshift uncertainty. The supplied master
retains its adopted redshift hierarchy.

The two cited reference CSVs in [data/literature](data/literature/README.md) support
the automatic literature crossmatch. Confirmed dual comparisons retain
`0 <= z <= 4`, separation at most 1.5 arcsec and a 1-arcsec positional finder.
Removal requires verified pair identity and a confirmed contaminant
classification. Uncertain or conflicting classifications remain flagged.
Visual inspection is not read or applied.

## Chance-alignment controls

Copy `chance_alignment.example.json` to `chance_alignment.local.json` and set its
full quasar master and published Quaia reference paths. Run the complete sequence,
including the independent control experiment, with:

```sh
python -m scripts.run_analysis --chance-config chance_alignment.local.json
```

The chance-alignment calculations follow candidate vetting and use their existing
public control downloads. They estimate method false-positive rates from the
input populations, not contamination within the literature-cleaned final list.
Their separate workspace contains the results. To run them after candidate vetting:

```sh
python -m scripts.run_chance_alignments --config chance_alignment.local.json --resume
```

## Verification

```sh
python -m unittest discover -s verification -p 'test_*.py'
python verification/verify_redshift_uncertainty.py
python verification/verify_input_contracts.py --models-only
```

Tests use local fixtures and mocked acquisition; they do not connect to WSDB or
perform a full survey run. After running the analysis, inspect its consistency:

```sh
python verification/verify_example_run.py
```

Automated tests establish implementation behavior. A complete run on the supplied
survey data is required to verify the actual candidate count. A replacement
master or survey release can change discovery membership, calibration and counts.
