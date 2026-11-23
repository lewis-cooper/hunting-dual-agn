# Chance-alignment controls

This package tests the acceptance of unrelated source pairs using the
production QSO scorers. The modules in `engine/` implement catalogue queries,
matching, scoring and statistical summaries. `prepare_inputs.py` prepares the
control inputs, and `scripts/run_chance_alignments.py` and
`scripts/run_analysis.py` run the calculations in sequence. Generated controls
and outputs are kept in a separate ignored workspace.

The measured quantity is **P(pass | unrelated matched control, known outcome)**,
for the brightness/redshift/sky distribution of the pre-vetting inputs.
It is not a selected-sample contaminant count or purity. The literature
crossmatch is maintained separately and is not included in this repository.

## Complete sequence

Install the dependencies and configure the input/model paths in the
[main README](../README.md). Copy `chance_alignment.example.json` to
`chance_alignment.local.json` and supply the full quasar master and Quaia paths.

```sh
python -m scripts.run_analysis --chance-config chance_alignment.local.json --resume
```

This runs the existing vetting stages and checks first, then the control stages.
The existing `run_frozen_analysis.sh` remains a vetting-only entry point, now
delegating to the same Python sequence so its stage order is maintained in one place.

## Control stages

```sh
python -m scripts.run_chance_alignments --config chance_alignment.local.json --stage prepare
python -m scripts.run_chance_alignments --config chance_alignment.local.json --stage fetch --resume
python -m scripts.run_chance_alignments --config chance_alignment.local.json --stage components --resume
python -m scripts.run_chance_alignments --config chance_alignment.local.json --stage blends --resume
python -m scripts.run_chance_alignments --config chance_alignment.local.json --stage report
```

Omit `--stage` to run all stages. `--dry-run` shows the connection without querying
or scoring. Paths are explicit through the config and the existing vetting
environment variables. Optional config keys are `input_directory`,
`products_directory`, `vetting_work_directory`, `model_directory` and
`dustmaps_directory`; their defaults come from `settings.py`. The required
master/Quaia schema is in the [input guide](../data/README.md).

1. **Prepare:** derive pre-vetting input descriptors and empirical DCR/morphology
   factors from the repo's own vetting outputs. Sample nested 30,000 and 60,000
   real quasar fields, with the original cell allocation, weights and seed.
2. **Fetch:** query Gaia DR3 and Legacy DR9 to 60 arcsec, associate sources and
   identify isolated quasar donors. Completed 30,000-field cones are reused in
   the expanded sample after coordinate checks.
3. **Components:** match actual input primaries to measured field sources and
   apply the production component scorer. Keep the sufficient 30,000-field
   one-Gaia run; resolved and two-Gaia imaging use the expanded library.
4. **Blends:** obtain independent native fluxes, use isolated matched quasar
   donors, and score their flux sums. Preserve the original 100 Legacy / 20
   fallback mock allocation, targeted coverage extension and Legacy expansion.
5. **Report:** assemble batches, retain shared source/field identities, calculate
   weighted rates and 95% uncertainty, and write tables, figures and `RESULTS.md`.

The maximum library has 60,000 real catalogue fields. The test pairs are
constructed associations between real sources, not 60,000 observed close pairs.
Flux sums do not constitute image injection or close-pair profile fitting.

## Outputs and resuming

The default workspace is `chance_alignment_work/`, ignored by Git. Prepared
input descriptors and calibration provenance are under `control_inputs/`.
Field/query/score checkpoints and logs are under `run/`. Combined outputs are
under `run/paired_60000/`:

- `repaired_summary_by_route.csv` and `repaired_summary_by_latitude.csv`;
- individual component/blended score Parquets and matching audits;
- magnitude-matching and separate outer-neighbour density diagnostics;
- `repaired_control_pass_rates.pdf` / `.png` and `RESULTS.md`.

Files with `repaired` in their names contain the matched-control results.
No additional selected-candidate lists are written. The six existing vetting
products are unaffected. `--resume` reuses completed steps and cached queries;
changed input contents or settings require a new `workspace_directory`.

## Offline replay and verification

With the separately supplied control-data archive unpacked:

```sh
python -m scripts.run_chance_alignments --replay-data /path/to/control_data --verify-reference
python -m unittest discover -s verification -p 'test_*.py'
```

The replay recalculates the statistics from individual control rows, checks their
weights and identities, and compares five numerical result tables to the supplied
reference. It needs no production model assets or survey queries. A fresh rebuild
does not need this archive: its input descriptors are generated from vetting and
its field lists are sampled from the supplied full catalogues.

Standalone utilities are available under this package: `run_controls.py`,
`reproduce_results.py` and `verify_models.py`. They accept explicit configuration
and control-data paths. `verify_models.py` uses the small
real scoring examples in the optional control-data archive and the corresponding
standalone config; it is not needed to run the full repo sequence.

## Preserved scientific rules and qualifications

- The threshold is read from the repo's calibration product, currently
  ln R >= -6.81. Either model may pass; imaging components must pass under the
  same model variant. Blends require finite QQ/QS scores in both and a positive
  score in either.
- Secondary PM must be below 3 when finite; missing PM is retained and flagged.
  Imaging primary PM must be finite and below 3. Quaia uncertainties use the
  existing production likelihood treatment.
- DCR applies only to one-Gaia imaging and is represented by a survival fraction
  from the real cut flow. It is not directly measured on the constructed random
  pair geometries. This transfer assumption is unchanged.
- Native fluxes/variances are added independently; missing fluxes are never zero
  filled. Unknown outcomes appear in separate sensitivity ranges, not as failures.
- Reused sources and donor fields remain linked in the uncertainty calculation.
  Raw mock combinations are not independent measurements.
- The wider 10-30 arcsec annulus is a separate check. Wide photometry, sky matching,
  available band masks and empirical pair factors remain transfer assumptions.
- The 10-60 arcsec density tables are separate geometric diagnostics. They do not
  produce a selected-sample contamination estimate or distinguish lenses from duals.

No personal credentials, WSDB code, input catalogues, control archives or generated
figures are included in Git.
