# Supplied local observations

The candidate workflow reads local observations and computes its own discovery
and measurement tables. Data acquisition is the user's responsibility. It does
not connect to a survey database, download catalogues or fetch images.

## Quasar and Gaia discovery inputs

`configs/quasar_catalogue.example.json` maps the supplied quasar catalogue.
It needs a unique object ID, RA/Dec in degrees and redshift. Supply either a
catalogue-wide redshift type or per-row types. Photometric redshifts require
finite nonnegative errors; spectroscopic values receive zero additional error.
Mixed-catalogue labels and JSON native-ID/error fields are supported.

`gaia_catalogue` must contain the complete relevant Gaia detections with exact
`source_id`, `ra`, `dec` and `phot_g_mean_mag`. It must include non-quasar
companions and isolated fields, rather than only known pairs. Supply
`gaia_coverage` with one row per input quasar: `object_id`, `ra_deg`, `dec_deg`,
`query_radius_arcsec >= 2`. These records certify the supplied observations;
an absent source is meaningful only within complete coverage.

The adapter matches all Gaia detections within 1 arcsec of the quasars and counts
Gaia sources within 1 arcsec of each possible Gaia anchor. Counts include the
central source. Discovery retains `cnt1 > 1`, chooses the nearest qualifying Gaia
match per quasar, then applies the historical catalogue-position deduplication.
Final physical-pair deduplication remains in the selection/combination stages.
Neither input rows nor candidate counts are required to equal 712.

### Historical survey selection

To reproduce the original Gaia-neighbour input definition, supply associations
from SDSS DR16Q (`is_qso_final = 1`), the original DESI DR1 quasar catalogue
(`vasily.desi_dr1_qso_zpix_v1` in the historical database), and Quaia G20 v1.
The original DESI query imposed no additional `ZCAT_PRIMARY`, `ZWARN` or
`SPECTYPE` restriction on that quasar table. A differently selected DESI
catalogue is a different input population.

Keep **all** SDSS/DESI Gaia DR3 associations within 1 arcsec, and use Quaia's
published Gaia IDs directly. Preserve the survey quasar coordinates in `qso_ra`
and `qso_dec`; substituting Gaia positions can change duplicate grouping.
Supply the associations and EDR3 auxiliary neighbour counts through
`configs/one_arcsecond_search.example.json`. Selection applies `cnt1 > 1`
before choosing the nearest qualifying match. Groups use a 1-arcsec positional
link and the original representative priority SDSS, DESI, then Quaia. This
grouping priority is separate from later redshift adoption.

The workflow does not add a `ZCAT_PRIMARY` requirement or require a particular
historical count. The supplied input definition and original identities must be
compared to establish reproduction; agreement in total counts alone is
insufficient. The combined imaging master does not define this historical
Gaia-neighbour input population.

## Imaging discovery inputs

`configs/imaging_search.example.json` describes the supplied Gaia-linked
master and local survey extracts. Preserve master `qso_parent_index` values:
Gaia/CCD query records must refer to that same population.

| Input | Required fields |
| --- | --- |
| Master | `qso_parent_index`, exact Gaia `source_id`, `gaia_ra`, `gaia_dec`, `adopted_redshift`, `adopted_redshift_catalogue`; redshift error/type when supplied |
| Tractor | `release`, `brickid`, `objid`, `ra`, `dec`, `type`, `maskbits`, `fracflux_r`; g/r/z fluxes, inverse variances and observation counts |
| Gaia | Source ID, RA/Dec, G magnitude, PM components/errors/correlation, `parallax_over_error`, `ruwe` |
| g-band CCDs | Centre, camera/filter/quality, dimensions, TAN WCS, observation MJD, exposure time and `sig1` |

Tractor coverage records are keyed by master Gaia `source_id`, with `ra_deg`,
`dec_deg` and `query_radius_arcsec >= 10`. Gaia component coverage contains cone
centres and radii of at least 0.7 arcsec. CCD coverage uses `qso_parent_index`,
centre and `query_radius_deg >= 0.65`. Records must include successful empty
fields. Matching requires the coverage coordinates to agree with the data.

The discovery searches the supplied Tractor sources within 1.5 arcsec, validates
two positive-flux PSF components, associates Gaia detections and reconstructs DCR
axes from local exposure metadata. North/south resolution, colour-based primary
assignment and survey-selection rules follow the imaging baseline introduced in
commit `509154f`.

The generated channel, component, labelled, Tractor and Quaia-reference FITS
tables are placed under `PAPER1_WORK_ROOT/imaging_search/inputs`.

## Component measurements and images

`configs/measurement_inputs.example.json` supplies Gaia/Tractor data for
the discovered Gaia-neighbour systems, exact DESI fallback photometry where
needed, PS1/NSC/SkyMapper photometry, and the image-path template. Native column
maps are canonical-to-file mappings. Fluxes are nanomaggies, inverse variances
are inverse squared nanomaggies, and positions are ICRS degrees.

Positional survey extracts need coverage keyed by exact Gaia anchor `source_id`
with centre and radius. Set `full_catalogue: true` only for a complete catalogue
covering the relevant region. It is not a way to mark a partial extraction as
complete. A supplementary survey may be explicitly disabled; that removes its
bands from scoring and may change the result.

PS1 stack catalogues can contain several measurements with the same object ID.
For those extracts, set `record_identity` to the unique native-row fields, for
example `["source_id", "ra", "dec"]` when that combination is unique. The
builder selects the nearest native row within 1 arcsec, matching the original
PS1 query, and records its coordinates. It does not keep an arbitrary first
object-ID duplicate. Conflicting rows with the same declared identity stop the
run. The default identity remains `source_id` for single-row catalogues.

The builder calculates component association, PSF/SNR/quality flags, structural
resolution classes and per-band eligibility. Shared nearest Tractor detections
are withheld from resolved component photometry.

`morphology_images.path_template` locates supplied FITS images. It may use row
fields such as `qso_parent_index`. For multiband cubes specify `r_band_index`.
The measurement uses image WCS and component coordinates, baseline-subtracted
parallel/perpendicular profile moments, 1-arcsec strip half-width, 1.5-arcsec
margins and one-pixel bins. It calculates the width ratio before applying the
`R_perp >= 1.05` cut. Missing required images stop the workflow.

## Scope and reproducibility

The package accepts supplied files and applies the analysis; obtaining those
files is a separate task. A quasar coordinate/redshift list alone does not
contain the observations needed for either complete selection pipeline.

The supplied population, releases, coverage and rounding conventions determine
discovery membership and calibration. The historical auxiliary neighbour-count
table and a newly supplied Gaia extract must be compared rather than assumed
identical. The code does not enforce historical candidate counts.

Model bundles and dust maps are separately installed local prerequisites.
Visual inspection is outside the workflow. The two cited literature CSVs are
included; [data/literature/README.md](data/literature/README.md) explains their matching and
classification policy. The optional chance-alignment experiment is separate
from the local candidate workflow and retains its own public control acquisition.

## Publication catalogue

The automatic catalogue in `products/final_candidates.csv` is distinct from the
curated publication files in `catalogue/`. `data/publication/exclusions.csv`
records the externally identified corrupted-data artefact, with exact Gaia
identities and its exclusion reason. `export_catalogue.py` applies that explicit
decision after selection, preserves original names and stable IDs, and writes
the combined publication catalogue and its two tier tables. No scoring or
selection threshold is changed. The Euclid visual-review inputs are described
in `data/publication/README.md`; image cutouts and MOC FITS files remain supplied
local observations.
