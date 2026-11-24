# Separately supplied analysis inputs

Survey catalogues, visual-inspection decisions and image data are not included
in this release. The one-arcsecond route operates on prepared tables. The
imaging route finds its initial pairs from a supplied Gaia-linked quasar master
and local survey measurements. Set `PAPER1_INPUT_ROOT` to a directory containing
the one-arcsecond and other external vetting inputs below.
This directory is ignored by Git when placed at `data/inputs/`.

The two small, cited literature reference snapshots in
[`literature/`](literature/README.md) are included in the repository and are
separate from these external survey inputs.

## Catalogue and photometry files

| Files | Role and keys |
| --- | --- |
| `original_parent_712.csv`, `eligibility_audit_712.csv`, `secondary_measurements_324.csv` | The 712 parent systems, eligibility measurements and 324-entry original resolved subset. Join on exact `qso_group_id`. |
| `one_arcsecond_resolution.csv`, `one_arcsecond_components.fits`, `one_arcsecond_tractor.csv`, `desi_target_photometry.csv` | Resolution assignments, component morphology, Tractor measurements and DESI fallback photometry. Tractor identity is `release`, `brickid`, `objid`. |
| `quaia_712_reference.fits` | Published Quaia `source_id`, `redshift_quaia` and `redshift_quaia_err` for the adopted one-arcsecond parent redshifts. |
| `ps1_measurements.csv`, `nsc_measurements.csv`, `skymapper_measurements.csv` | Blended-system photometry with native magnitudes, errors and quality flags. Join by `qso_group_id`. |
| `sdss_measurements.csv`, `allwise_measurements.csv`, `vhs_measurements.csv` | Inputs to the photometry inventory stage. |

Positions are ICRS degrees; separations are arcseconds. Tractor fluxes are
nanomaggies and inverse variances are inverse squared nanomaggies. Supply
native photometry: the scoring code applies extinction correction and converts
to luptitudes. Preserve Gaia and catalogue identifiers as exact strings or
integers. Redshift errors are dimensionless.

These are prepared associations with the columns used by the preparation and
scoring scripts for the one-arcsecond route. Obtain survey measurements from their public archives,
and perform the required joins and quality handling before running this code:

- [SDSS DR16Q](https://www.sdss4.org/dr16/algorithms/qso_catalog/)
- [DESI DR1](https://data.desi.lbl.gov/doc/releases/dr1/)
- [Quaia](https://irsa.ipac.caltech.edu/data/Quaia/overview.html)
- [Gaia Archive](https://gea.esac.esa.int/archive/)
- [Legacy Surveys DR9](https://www.legacysurvey.org/dr9/)
- [Pan-STARRS](https://panstarrs.stsci.edu/), [NSC](https://datalab.noirlab.edu/dataexploration/dataaccess/nsc) and [SkyMapper](https://skymapper.anu.edu.au/data-release/dr4/)

## Imaging candidate-search inputs

Set local paths in `imaging_search.local.json`, copied from the root example.
All files are local CSV or FITS catalogues. Tractor inputs may be a list of
files/globs; FITS files are read in chunks and CSV files are streamed. No
catalogue acquisition or network query is performed by this stage.

| Config key | Required measurements |
| --- | --- |
| `master_catalogue` | Gaia-linked master with unique `qso_parent_index`, exact `source_id`, `gaia_ra`, `gaia_dec`, `adopted_redshift`, `adopted_redshift_catalogue`. Optional `sdss_name`, `q_ra`, `q_dec`, catalogue memberships and original provenance are retained. Missing query coordinates default to the Gaia coordinates. |
| `tractor_catalogues` | Local DR9 sources with `release`, `brickid`, `objid`, `type`, `ra`, `dec`, `maskbits`, `fracflux_r`, and `flux_*`, `flux_ivar_*`, `nobs_*` for g/r/z. `brick_primary` is used when present. Raw unassociated source catalogues are spatially matched locally. Existing field-neighbour caches are also supported when they include exact `parent_gaia_source_id`; old numerical parent indices alone are insufficient. |
| `tractor_coverage` | One successful-query record per master Gaia `source_id`, with `ra_deg`, `dec_deg`, `query_radius_arcsec` (at least 10). Include successful empty queries. The coordinates must agree with the master. |
| `gaia_catalogue` | Local Gaia DR3 sources with exact `source_id`, `ra`, `dec`, `phot_g_mean_mag`, `pmra`, `pmdec`, `pmra_error`, `pmdec_error`, `pmra_pmdec_corr`, `parallax_over_error`, `ruwe`. Repeated source IDs from overlapping complete cones are deduplicated. |
| `gaia_coverage` | Successful component-cone centres `ra_deg`, `dec_deg`, `query_radius_arcsec` (at least 0.7), including cones with no detections. |
| `ccd_catalogue` | DR9 g-band CCD metadata: `camera`, `filter`, `ccd_cuts`, `mjd_obs`, `exptime`, `sig1`, `width`, `height`, `crpix1/2`, `crval1/2`, `cd1_1`, `cd1_2`, `cd2_1`, `cd2_2`. Unassociated CCD catalogues also need centre `ra`, `dec`; field-query caches instead use `qso_parent_index` and may retain `target_ra`, `target_dec` for coordinate checks. |
| `ccd_coverage` | Successful query records `qso_parent_index`, `ra_deg`, `dec_deg`, `query_radius_deg` (at least 0.65), centred on the discovered r-band primary. Include successful empty queries. |
| `quaia_catalogue` | Published `source_id`, `redshift_quaia`, `redshift_quaia_err`. A complete superset is acceptable; only the discovered master records are used, and adopted Quaia redshifts must agree exactly within 1e-12. |

Coverage records are data-acquisition logs, not candidate selections. Their
purpose is to distinguish a genuinely empty cone from an unqueried field.
They must describe complete supplied catalogue extracts; the presence of one
nearby source does not establish coverage. For bulk survey files, establish
which complete cones their footprints cover before recording those fields as
available. Missing or inconsistent coverage stops the search rather than
treating those fields as non-detections.

The initial Tractor search uses a 1.5 arcsec aperture and retains the original
full-master nearest-source colour reference. Do not apply a new latitude or
magnitude cut to this master if reproducing the paper's counts. Later vetting
applies the latitude, photometric, PM and DCR requirements. Component roles
are fixed using the original local g-r median within z +/- 0.1, excluding the
target and requiring at least 20 reference objects; otherwise the nearest
component to the supplied quasar position is used, with the original r-band
fallback retained when necessary.

Master CSVs use the same standard pandas floating-point parsing as the original
search. This preserves the adopted redshift values exactly: the model's
support-sampling seed hashes its numerical inputs, so even a last-bit change
from a different CSV parser can change a borderline support decision. Other
catalogue CSVs retain round-trip parsing. This convention is checked when
comparing the rebuilt inputs and scores with the reference run.

With only the master/Tractor inputs available, `--stage discover` writes the
Gaia and CCD query-position tables for subsequent data acquisition. Supply
the additional local files later and rerun the complete stage. Search and
association outputs are generated under `PAPER1_WORK_ROOT/imaging_search/`;
the five vetting inputs are in its `inputs/` subdirectory:

| Generated file | Role |
| --- | --- |
| `imaging_channel.fits` | Initial ranked Tractor pairs, including the PSF–PSF discovery sample. |
| `imaging_components.fits` | Validated, r-ranked PSF components and exact Tractor identities. |
| `imaging_labelled.fits` | Fixed primary/secondary roles, component Gaia associations, covariance-aware PM significance and exposure-based DCR angle. |
| `imaging_tractor.fits` | Local component photometry and quality metadata. |
| `quaia_imaging_reference.fits` | Published Quaia redshifts/errors for the discovered master records. |

Separate Gaia detections are not forced: if the same Gaia source is the nearest
match to both components, the system belongs to the one-distinct-Gaia tier.
Missing Gaia PM remains missing for the following selection policy. Images,
image-model decomposition and morphology-profile generation are outside this
search stage; the supplied morphology table below remains required.

## Imaging morphology measurements

Supply a CSV with one row per measured `qso_parent_index` and these columns:

| Column | Meaning |
| --- | --- |
| `qso_parent_index` | Exact imaging parent identifier; unique within the table. |
| `image_axis_sigma_arcsec` | Previously measured along-pair image width in arcseconds. |
| `image_perp_sigma_arcsec` | Previously measured perpendicular image width in arcseconds. |
| `r_perp_sigma_ratio` | Previously measured width ratio used for the $R_\perp\geq1.05$ selection. |
| `morphology_measurement_source` | Provenance label for the supplied measurement. |

The default location is `PAPER1_INPUT_ROOT/imaging_morphology_measurements.csv`.
Override it with `PAPER1_IMAGING_MORPHOLOGY_PATH`.
The selector joins measurements by parent identifier and rejects candidate
systems without a finite ratio. Image-profile generation and its input images
are outside this release.

## Visual-inspection decisions

Supply `visual_inspection_rejections.csv` with columns
`one_arcsecond_qso_group_id`, `imaging_qso_parent_index` and `reason`.
Use an empty identifier field when a rejection applies only to the other
search. The file records the decisions made for the sample being analysed;
it is not a substitute for inspecting new candidates.

The default location is under `PAPER1_INPUT_ROOT`. Override it with
`PAPER1_VISUAL_REJECTIONS_PATH`. The paper's original decisions are needed to
reproduce its reference counts.

## External dependencies

Install the pinned model package and SFD dust maps using the
[main README](../README.md). `verification/verify_input_contracts.py` checks
that the 16 external vetting inputs and five generated imaging tables are present, the model manifests/revision match,
and the extinction maps exist. File availability does not establish that a
replacement catalogue reproduces the paper's scientific associations.

## Additional inputs for chance-alignment controls

The fresh control workflow starts from the vetting outputs and working tables;
these are generated by the sequence in the main README. In particular it reads
`one_arcsecond_scores.csv`, `imaging_scores.csv`, `threshold_calibration.json`
and the working `imaging/labelled_systems_audit.csv`, together with the supplied
generated component/Tractor FITS files above. It does not read a final candidate list.

Two additional, separately supplied catalogues are required:

| Config key | Required columns and meaning |
| --- | --- |
| `master_catalogue` | The Gaia-linked SDSS+DESI+Quaia master CSV: `qso_parent_index`, `source_id`, `gaia_ra`, `gaia_dec`, `phot_g_mean_mag`, `adopted_redshift`, `adopted_redshift_catalogue`. Positions are ICRS degrees. Preserve exact Gaia identifiers. |
| `quaia_catalogue` | The full published Quaia G<20 FITS catalogue: `source_id`, `redshift_quaia`, `redshift_quaia_err`. The reduced vetting reference files do not cover every possible control donor. |

The combined master is a prepared survey association table supplied by the
research analysis, rather than a single public catalogue download. Obtain the
SDSS, DESI, Quaia and Gaia inputs from their archives listed above and preserve
the adopted source associations/redshift hierarchy, or supply the original
prepared master to reproduce the reference experiment. The reference master
contains 1,157,489 entries at |b| >= 25 degrees; G<20 applies to Quaia, not to
all SDSS/DESI primaries.

The control scripts fetch their Gaia DR3, Legacy DR9 and fallback PS1/NSC/SkyMapper
measurements from public services, with successful queries cached locally.
They use no WSDB connection or personal credentials. Input catalogues, downloaded
controls and generated outputs stay outside Git. Paths in
`chance_alignment.local.json` are resolved relative to that config file.
