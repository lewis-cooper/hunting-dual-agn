# Chance-alignment controls

This package tests the acceptance of unrelated source pairs using the
production QSO scorers. The modules in `engine/` implement catalogue queries,
matching, scoring and statistical summaries. `prepare_inputs.py` prepares the
control inputs, and `run_chance_alignments.py` and
`run_analysis.py` run the calculations in sequence. Generated controls
and outputs are kept in a separate ignored workspace.

The measured quantity is **P(pass | unrelated matched control, known outcome)**,
for the brightness/redshift/sky distribution of the pre-vetting inputs.
It is not a selected-sample contaminant count or purity. The separate
literature stage is documented in the [main README](../README.md); it removes
verified, confirmed contaminants from the candidate catalogue without changing
the populations or scorers used by these method-level controls.

## Complete sequence

Install the dependencies and configure the input/model paths in the
[main README](../README.md). Copy `configs/chance_alignment.example.json` to
`chance_alignment.local.json` and supply the full quasar master and Quaia paths.

```sh
python run_analysis.py --stage all --pipeline both --chance-config chance_alignment.local.json --resume
```

This runs local candidate discovery, vetting/literature stages and
checks first, then the control stages. The control adapter reads the generated
imaging component/Tractor tables from the search workspace. The control
scorers and statistical calculations are unchanged.

## Control stages

```sh
python run_chance_alignments.py --config chance_alignment.local.json --stage prepare
python run_chance_alignments.py --config chance_alignment.local.json --stage fetch --resume
python run_chance_alignments.py --config chance_alignment.local.json --stage components --resume
python run_chance_alignments.py --config chance_alignment.local.json --stage blends --resume
python run_chance_alignments.py --config chance_alignment.local.json --stage report
```

Omit `--stage` to run all stages. `--dry-run` shows the connection without querying
or scoring. Paths are explicit through the config and the existing vetting
environment variables. Optional config keys are `input_directory`,
`products_directory`, `vetting_work_directory`, `model_directory` and
`dustmaps_directory`; their defaults come from `dual_agn/settings.py`. The required
master/Quaia schema is in the [input guide](../DATA.md).

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
   weighted rates and 95% uncertainty, compare detected close neighbours with
   the 10–60 arcsec density in the same fields, and write tables, figures and
   `RESULTS.md`.

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
- `matched_density_comparison.csv`: observed close counts and geometric random
  predictions at several radii, including checks using 10–30 and 30–60 arcsec;
- `matched_density_field_counts.csv`: per-field counts, exact Gaia identities
  and population weights, including fields with no neighbours;
- `observed_vs_field_density.pdf` / `.png` and `FIELD_DENSITY.md`;
- `field_density_inputs/`: the source catalogue, central Gaia IDs, complete
  close/outer Gaia detections, coverage and query metadata needed for replay;
- `repaired_control_pass_rates.pdf` / `.png` and `RESULTS.md`.

Files with `repaired` in their names contain the matched-control results.
No additional selected-candidate lists are written. The six existing vetting
products are unaffected. `--resume` reuses completed steps and cached queries;
changed input contents or settings require a new `workspace_directory`.

## Offline replay and verification

With the separately supplied control-data archive unpacked:

```sh
python run_chance_alignments.py --replay-data /path/to/control_data --verify-reference
python -m unittest discover -s tests -p 'test_*.py'
```

The replay recalculates the statistics from individual control rows, checks their
weights and identities, and compares five numerical result tables to the supplied
reference. It needs no production model assets or survey queries. A fresh rebuild
does not need this archive: its input descriptors are generated from vetting and
its field lists are sampled from the supplied full catalogues.

New control-data bundles need a `field_density/` directory containing the
four files exported under `run/paired_60000/field_density_inputs/`:
`master_targets.csv`, `gaia_neighbours.parquet`, `cone_coverage.csv` and
`query_configuration.json`. The Gaia table includes the central sources so
they can be excluded by exact identity, rather than an assumed position cut.
Keep this data directory separate from the replay output.

Older control archives omit the central Gaia IDs and coverage metadata.
Supply the completed source cones explicitly when replaying those archives:

```sh
python run_chance_alignments.py --replay-data /path/to/control_data \
  --density-source /path/to/run/fields_60000 --verify-reference
```

The replay exports the density inputs again, regenerates the matched comparison
and includes it in `RESULTS.md`. The five historical reference tables are checked
as before; the new comparison additionally validates source identities and cone
coverage. For a density-only calculation from cached cones:

```sh
python -m chance_alignments.engine.analyze_field_density \
  --source /path/to/run/fields_60000 --output /path/to/density_results \
  --export-inputs /path/to/control_data/field_density
```

The public control entry point is `run_chance_alignments.py`.
`run_controls.py` implements its execution sequence and `reproduce_results.py`
implements offline replay.

## Preserved scientific rules and qualifications

- The threshold is read from the current run's calibration product.
  Either model may pass; imaging components must pass under the
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
- The matched comparison predicts close counts by scaling outer counts by
  `theta**2 / (60**2 - 10**2)`, assuming complete detection. Both the raw sample
  counts and the inclusion-weighted population estimates are reported. The
  observed/predicted ratio measures neither the unrelated fraction nor detection
  completeness: missed close detections and associated companions are not
  independently constrained. Agreement alone would not show that most input pairs
  are chance alignments. Do not combine pooled densities and pass rates into a
  claimed actual contaminant count. Bootstrap intervals exclude completeness and
  population-transfer uncertainty.

No personal credentials, input catalogues, control archives or generated
figures are included in Git.
