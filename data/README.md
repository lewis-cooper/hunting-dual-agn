# Separately supplied analysis inputs

Input catalogues, visual-inspection decisions and image data are not included
in this release. The scripts operate on prepared tables from the research
analysis. Set `PAPER1_INPUT_ROOT` to a directory containing the files below.
This directory is ignored by Git when placed at `data/inputs/`.

## Catalogue and photometry files

| Files | Role and keys |
| --- | --- |
| `original_parent_712.csv`, `eligibility_audit_712.csv`, `secondary_measurements_324.csv` | The 712 parent systems, eligibility measurements and 324-entry original resolved subset. Join on exact `qso_group_id`. |
| `one_arcsecond_resolution.csv`, `one_arcsecond_components.fits`, `one_arcsecond_tractor.csv`, `desi_target_photometry.csv` | Resolution assignments, component morphology, Tractor measurements and DESI fallback photometry. Tractor identity is `release`, `brickid`, `objid`. |
| `quaia_712_reference.fits` | Published Quaia `source_id`, `redshift_quaia` and `redshift_quaia_err` for the adopted one-arcsecond parent redshifts. |
| `ps1_measurements.csv`, `nsc_measurements.csv`, `skymapper_measurements.csv` | Blended-system photometry with native magnitudes, errors and quality flags. Join by `qso_group_id`. |
| `sdss_measurements.csv`, `allwise_measurements.csv`, `vhs_measurements.csv` | Inputs to the photometry inventory stage. |
| `imaging_channel.fits`, `imaging_components.fits`, `imaging_labelled.fits` | Imaging-search pairs, fixed component assignments, Gaia associations, PM significance and folded pair-DCR angle. Join on `qso_parent_index`. |
| `imaging_tractor.fits` | Component photometry matched by exact Tractor identity. |
| `quaia_imaging_reference.fits` | Published Quaia redshifts and uncertainties, matched on the master record's `source_id` supplying the adopted redshift. |

Positions are ICRS degrees; separations are arcseconds. Tractor fluxes are
nanomaggies and inverse variances are inverse squared nanomaggies. Supply
native photometry: the scoring code applies extinction correction and converts
to luptitudes. Preserve Gaia and catalogue identifiers as exact strings or
integers. Redshift errors are dimensionless.

These are prepared associations with the columns used by the preparation and
scoring scripts. Obtain the survey measurements from their public archives,
and perform the required joins and quality handling before running this code:

- [SDSS DR16Q](https://www.sdss4.org/dr16/algorithms/qso_catalog/)
- [DESI DR1](https://data.desi.lbl.gov/doc/releases/dr1/)
- [Quaia](https://irsa.ipac.caltech.edu/data/Quaia/overview.html)
- [Gaia Archive](https://gea.esac.esa.int/archive/)
- [Legacy Surveys DR9](https://www.legacysurvey.org/dr9/)
- [Pan-STARRS](https://panstarrs.stsci.edu/), [NSC](https://datalab.noirlab.edu/dataexploration/dataaccess/nsc) and [SkyMapper](https://skymapper.anu.edu.au/data-release/dr4/)

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
that the 21 local input files are present, the model manifests/revision match,
and the extinction maps exist. File availability does not establish that a
replacement catalogue reproduces the paper's scientific associations.
