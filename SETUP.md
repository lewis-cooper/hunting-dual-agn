# Setup and execution

## Install


Use Python 3.11 or newer. From this directory:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
git clone https://github.com/vasilybelokurov/qso_p_color.git ../qso_p_color
git -C ../qso_p_color checkout dbc392623d27e22744d5224f9f1a3fd52f6230a8
python -m pip install -e ../qso_p_color
export QSO_PCOLOR_PACKAGE_ROOT="$(cd ../qso_p_color && pwd)"
```

The pinned package supplies both model bundles in `models/multisurvey_psf/`.
Obtain the two [SFD dust maps](https://doi.org/10.7910/DVN/EWCNL5), place them
under `sfd/`, and set `QSO_PCOLOR_DUSTMAPS_DIR` to its parent directory.
The workflow checks these local prerequisites before scoring.

## Supply the observations


Copy the templates into the repository root:

```sh
cp configs/quasar_catalogue.example.json quasar_catalogue.local.json
cp configs/imaging_search.example.json imaging_search.local.json
cp configs/measurement_inputs.example.json measurement_inputs.local.json
```

Configure their file locations and column maps. Paths are resolved relative to
each configuration file. CSV and FITS catalogues are supported.

| Configuration | Supplied observations |
| --- | --- |
| `quasar_catalogue.local.json` | Quasar IDs, positions, redshifts/types/errors; complete Gaia neighbourhoods and coverage around the starting population |
| `imaging_search.local.json` | Gaia-linked imaging master, local Tractor and Gaia sources, g-band exposure metadata and coverage |
| `measurement_inputs.local.json` | Gaia/Tractor component measurements, available supplementary photometry, and local r-band FITS image locations |

The quasar column map supports other catalogues and sample sizes without Python
edits. Photometric redshifts require uncertainties. Native DESI IDs are needed
when requesting the corresponding target-photometry fallback.

The Gaia-neighbour route needs ordinary Gaia detections as well as quasars.
Complete 2-arcsec cones around the input quasars enclose a possible 1-arcsec-offset
Gaia match and its full 1-arcsec neighbourhood. The package computes counts;
users do not supply a selected parent sample or a neighbour-count summary in
this raw-Gaia mode.

Alternatively, supply quasar–Gaia associations and neighbour counts exported
from a survey database. The existing `cnt1` counts include the central source;
`cnt1 > 1` means at least one companion within 1 arcsec. Copy
`configs/one_arcsecond_search.example.json`, set its file paths, and run:

```sh
python run_analysis.py --pipeline both \
  --one-arcsecond-search-config one_arcsecond_search.local.json
```

This bypasses raw-Gaia neighbour counting. The workflow applies `cnt1 > 1`, chooses
the nearest qualifying association per quasar and deduplicates the quasar entries.
The exports must cover the desired input population. `saved_qualifying_matches`
accepts complete exports already filtered by `cnt1 > 1`; `survey_gaia_matches`
accepts unfiltered associations and corresponding counts. Neither mode reads a
previously grouped candidate list or requires a particular candidate count.
Native component photometry and FITS images are still needed for vetting.
Acquiring these exports is separate from this offline analysis package.

For the paper's original Gaia-neighbour selection, supply the original survey
associations, preserving survey coordinates and all Gaia matches. Do not replace
them with a combined master that changes DESI selection or duplicate-grouping
positions. No extra `ZCAT_PRIMARY=True` restriction is applied by discovery.
The exact historical input definition is described in
[DATA.md](DATA.md#historical-survey-selection).

The supplied imaging master and its parent IDs are preserved so that existing
Gaia/CCD coverage records retain their meaning. See [DATA.md](DATA.md)
for fields, units and coverage requirements.

## Run the candidate workflow


Check supplied files without searching or scoring:

```sh
python run_analysis.py --pipeline both --check-inputs
```

Run either pipeline or both:

```sh
python run_analysis.py --pipeline gaia_neighbour
python run_analysis.py --pipeline imaging
python run_analysis.py --pipeline both
```

`both` is the default. The Gaia-neighbour pipeline includes resolved and blended
routes. Imaging-only scoring automatically builds the Gaia-neighbour calibration
sample, but its final catalogue contains imaging candidates only. It does not
score the Gaia blended route.

Run only discovery to inspect initial counts before scoring:

```sh
python run_analysis.py --pipeline gaia_neighbour --stage discovery
python run_analysis.py --pipeline imaging --stage discovery
python run_analysis.py --pipeline both --stage discovery
```

Imaging-only discovery uses its supplied imaging master directly and needs no
Gaia-neighbour calibration data. The full scoring run requires that calibration
population. Calibration needs enough finite eligible resolved systems with
varying scores and PM significance; very small catalogues may not support it.

`--dry-run` prints the command sequence. `--resume` reuses discovery results only
when input/code/output hashes agree; scoring and selection run again. No candidate
count or fitted threshold is forced to a historical value.

Terminal output identifies each analysis step and reports elapsed time every
30 seconds while waiting for completion. Set `PAPER1_PROGRESS_INTERVAL_SECONDS`
to change that interval. A heartbeat indicates a live waiting process; it is not
a percentage-complete measurement or remaining-time estimate.

To keep test outputs separate:

```sh
export PAPER1_OUTPUT_ROOT=/path/to/test/products
export PAPER1_WORK_ROOT=/path/to/test/work
export PAPER1_FIGURE_ROOT=/path/to/test/figures
python run_analysis.py --pipeline both \
  --quasar-catalogue-config /path/to/test/quasar_catalogue.local.json \
  --imaging-search-config /path/to/test/imaging_search.local.json \
  --measurement-config /path/to/test/measurement_inputs.local.json
```

Completion prints the final candidate count and CSV location. Discovery reports
its own input and output counts in the working folder.

## Outputs


The combined run writes six files:

- `final_candidates.csv`
- `one_arcsecond_scores.csv`
- `imaging_scores.csv`
- `one_arcsecond_cutflow.csv`
- `imaging_cutflow.csv`
- `threshold_calibration.json`

A single-pipeline scoring run writes four: its score audit/cutflow, the threshold
and the final catalogue. Other intermediates stay under `PAPER1_WORK_ROOT` and
figures under `PAPER1_FIGURE_ROOT`. Use separate output folders for separate
pipeline runs. Identifiers should be read as strings.

The six products above and the contents of `figures/` are tracked by Git by
default. Working files and supplied inputs remain excluded. Analysis commands
write local outputs; publishing updated results uses the normal Git workflow:

```bash
git add products/ figures/
git commit -m "Update analysis results and figures"
git push
```

Commit products and figures together after the run has completed and its outputs
have been checked. The historical snapshot currently published is described in
[PROVENANCE.md](PROVENANCE.md).

## Optional method controls

See [chance_alignments/README.md](chance_alignments/README.md) for the separately requested control experiment.
