# Input data and generated measurements

The default workflow starts with a supplied quasar catalogue and acquires raw
survey records through the sources in `wsdb_inputs.local.json`. One configuration
contains the master column map, WSDB sources, the public exposure
sources and image-download settings. It builds all analysis-specific tables.

Positions are ICRS degrees. Query radii are arcseconds unless a name explicitly
ends in `_deg`. Tractor fluxes are nanomaggies and inverse variances are
nanomaggies^-2. Preserve identifiers as strings or native integer FITS columns;
never convert a large Gaia/DESI ID through floating point.

## Quasar catalogue

The inline `quasar_catalogue` object in `examples/wsdb_inputs.example.json` supports CSV
and FITS input and maps canonical fields to native columns:

| Field | Meaning |
| --- | --- |
| `object_id` | Unique, nonmissing input-row identity used for query coverage |
| `ra`, `dec` | Quasar coordinates in degrees |
| `redshift` | Adopted system redshift |
| `redshift_error` | Required for photometric redshifts; finite and nonnegative |
| `redshift_kind` | Optional per-row `spectroscopic` or `photometric` classification |
| `catalogue_label` | Optional per-row source-catalogue label |
| `name` | Optional display name |

Alternatively set scalar `catalogue_label` and `redshift_kind` in the config.
For a mixed master, `redshift_kind_by_catalogue` maps the original catalogue
labels to their types. `catalogue_labels` can map them to `SDSS_DR16Q`,
`DESI_DR1_QSO` and `QUAIA_G20_V1`, retaining the existing source-specific fallback
rules. Spectroscopic redshifts receive zero additional uncertainty.

`redshift_error_format: "json_singleton"` reads a single published error stored
in a JSON list for each photometric row. Spectroscopic rows may have empty error
lists. Ambiguous/missing photometric errors stop input validation.

A master may combine several native records. Map their JSON identifier fields
in `columns`, then choose them using `catalogue_object_id_columns`, keyed by the
original catalogue label, and `catalogue_object_id_format: "json_array"`. All
native IDs remain in the discovered parent. A single ID supplies the native
object identity; multiple IDs retain the master identity and explicit native-ID
list. DESI fallback photometry uses exact target IDs. Identical measurements from
several targets agree; conflicting target measurements are unavailable with a
`conflicting_native_targets` flag. No target is chosen arbitrarily. The original
single-target path retains its behavior.

The combined SDSS+DESI+Quaia master is a supported input rather than a bundled
catalogue. The supplied adopted redshift hierarchy is retained. Changing the
starting population or coordinates can change discovery membership; the workflow
does not force the historical 712 parent entries or a fixed final count.

## Data acquisition

`sources` maps a source name to `table`, `release` and `columns` (canonical to
native WSDB column names). The example maps the fields inspected in the WSDB
survey schemas. Only identifiers are accepted as table/column settings; query
expressions and credentials are rejected. Multiple tables can be supplied using
`{"tables": [{"table": "schema.north", "release": "north release", "columns": {...}},
{"table": "schema.south", "release": "south release", "columns": {...}}]}` with
identical ordered canonical fields. A source table must cover the intended sky
region and survey release. Query success alone does not establish survey footprint.

| Source | Canonical fields |
| --- | --- |
| `gaia` | `source_id`, `ra`, `dec`, `phot_g_mean_mag`, `pmra`, `pmdec`, their errors, `pmra_pmdec_corr`, `parallax_over_error`, `ruwe` |
| `tractor` | `release`, `brickid`, `objid`, `type`, `ra`, `dec`, `maskbits`, `fracflux_r`, plus g/r/z `flux_*`, `flux_ivar_*`, `nobs_*`; `brick_primary` is used when supplied |
| `desi` | `targetid`, `release`, g/r/z `flux_*` and `flux_ivar_*` |
| `ps1` | Source ID/position and native grizy Kron magnitudes/errors, mapped to `ps1_*kronmag`/`ps1_*kronmagerr` |
| `skymapper` | Source ID/position and uvgriz Petrosian magnitudes/errors, mapped to `skymapper_*_petro`/`skymapper_e_*_petro` |
| `nsc` | Source ID/position and ugrizyVR native magnitudes/errors, mapped to `nsc_*mag`/`nsc_*err` |
| `ccd` (optional WSDB alternative) | Centre position, `filter`, `ccd_cuts`, `camera`, dimensions, TAN WCS, `mjd_obs`, `exptime`, `sig1` |

The example uses `public_sources` for `ccd: "legacy_dr9_ccds"`. Data Lab
`ls_dr9.ccds_annotated` supplies quality-selected ground-based g exposures with
WCS and time metadata for DCR. Gaia epoch-photometry CCD tables cannot supply
these ground-based observations. The public inputs are recorded separately
from WSDB exports. They use no personal database credentials.

The acquisition sequence is:

1. Query complete 2-arcsec Gaia cones around every starting quasar. Include empty
   and isolated fields. Compute all 1-arcsec positional matches and cnt1 locally;
   cnt1 includes the central detection. Two-arcsec cones enclose a possible
   1-arcsec-offset anchor and its complete 1-arcsec neighbourhood.
2. Select cnt1 greater than one before choosing the nearest qualifying Gaia
   match per quasar, then deduplicate catalogue positions within 1 arcsec.
   Imaging separately uses the nearest Gaia anchor and one row per anchor.
3. Query Tractor sources to 10 arcsec around imaging anchors and any additional
   one-arcsecond anchors. The existing DR9 north/south selection and source
   identity rules are applied locally.
4. Discover PSF pairs within 1.5 arcsec. Query Gaia within 0.7 arcsec of each
   component and g-band CCD centres within 0.65 degrees of the r-ranked primary.
   CCD WCS then determines which exposures actually contain that position.
5. Query supplementary photometry in 2-arcsec anchor cones. The builders retain
   the original 1-arcsec matching aperture around the quasar coordinates. Query
   DESI fallback photometry by exact native target ID.
6. Build component tables, resolution classes and quality flags, then score and
   select using the existing rules. Measure r-band morphology from local images,
   downloading missing cutouts only when enabled.

Raw WSDB exports are batched in groups of at most 2,000 positions with individual
checkpoints. Query summaries record source release/mapping, positions/radii/code
fingerprints and output hashes. `--resume` reuses verified exports; it does not
check for later database updates. Overlapping batches can repeat native rows,
which the shared builders deduplicate. Public query responses have their own
hashed checkpoints. Only complete successful queries yield coverage records.
An error stops the workflow instead of becoming an empty field.

Generated configurations and raw files are under `PAPER1_WORK_ROOT/wsdb_inputs`.
Final products remain separate: six for both pipelines or four for a single
pipeline. The only server writes are temporary
query-position tables. Credential values are neither printed nor saved.

## Images and morphology

`morphology_images.download_missing: true` enables Legacy Surveys DR9 r-band
cutouts. The default is 100 pixels at 0.1 arcsec/pixel (10 arcsec), centred on
the supplied cutout/Gaia position with the original seven-decimal coordinate
rounding. Generic rows without either centre use the pair midpoint. These
settings match the original acquisition script. Downloaded FITS data must contain both component coordinates in their
WCS. Partial/invalid responses are not installed as images. Raw-image request
parameters and hashes are recorded alongside the image.

By default these files are cached under `PAPER1_WORK_ROOT/wsdb_inputs/rband_images`.
Set `path_template` to use another location; it can include `qso_parent_index`
or fields from the scored-system row. Set `r_band_index` for local multiband
cubes. Existing images are reused. Set `download_missing: false` for an offline
run; missing required images then stop the stage with example paths.

The existing morphology calculation forms baseline-subtracted profile moments
parallel/perpendicular to the pair at its midpoint, with 1-arcsec half-width,
1.5-arcsec margins and one-pixel bins. It writes the width ratio and image
provenance before the selection applies `R_perp >= 1.05`. No premeasured
morphology or visual-inspection table is read.

## Local-file route

Use `--data-source local` and the quasar, imaging and measurement configuration
examples. The master adapter needs raw Gaia detections and a coverage CSV keyed
by input `object_id`, with `ra_deg`, `dec_deg`, `query_radius_arcsec >= 2`.
It computes positional associations and counts itself.

The measurement builder reads native Gaia, Tractor, DESI, PS1, NSC and SkyMapper.
Their column maps are in `examples/local/measurement_inputs.example.json`. Positional inputs
need coverage keyed by exact anchor `source_id`, with cone centre and radius.
Successful empty cones must be included. Set `full_catalogue: true` only for a
complete survey catalogue covering the relevant region, never to bypass missing
queries. Shared nearest Tractor detections are withheld from resolved component
photometry, and per-band quality flags are calculated from native measurements.

Imaging Tractor coverage is centred on the master anchor with radius at least
10 arcsec. Gaia component coverage must enclose every 0.7-arcsec matching cone.
CCD coverage is keyed by `qso_parent_index` at the discovered r-ranked primary,
with centre and `query_radius_deg >= 0.65`. Changing discovery membership requires
rebuilding these coordinate-dependent associations.

Generated WSDB exports can be reused offline by passing their
`quasar_catalogue.json`, `imaging_search.json` and `measurement_inputs.json` to
the local runner. This recalculates the neighbour counts. Relocate absolute paths
when transferring these data. Raw extracts may be distributed separately where
survey terms permit; they are excluded from this Git repository.

## Models, controls and literature

Install the pinned model package and SFD maps as described in the main README.
`verification/verify_input_contracts.py --models-only` checks them without opening
survey data or a network connection. The full workflow also checks generated
input identities and output consistency. A useful calibration population and
adequate survey coverage are required for a replacement quasar catalogue.

Chance-alignment controls retain their separate public acquisition and scoring
sequence. Their config needs the full Gaia-linked quasar master and full published
Quaia catalogue, not just the reduced reference for selected systems. Controls
use input populations, not the final candidate list. The bundled literature
reference CSVs and their classification/matching policy are documented in
[literature/README.md](literature/README.md).
