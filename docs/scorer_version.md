# Scorer revision for the paper

This analysis requires the clean `qso_pcolor` checkout at
`dbc392623d27e22744d5224f9f1a3fd52f6230a8` (`dbc3926`). Model bundles being
unchanged does not make a later scorer interchangeable with this revision.

The paper uses the following behaviour:

- A source below the support cut (2/256), or outside both models, is rejected
  and receives no score. Later versions can score such sources and flag them.
- Legacy r below S/N 10 is rejected as `too_faint`. Later versions can select
  another reference band and retain the weak-r measurement in the likelihood.
- The workflow supplies the primary redshift. A later redshift-free
  `p_quasar` mode is outside this analysis.

Accepting previously rejected sources can change candidate membership and the
calibration population even if scores for previously accepted sources agree.
Using the later rules requires recomputing selection and threshold calibration,
then revising the method description. Euclid scoring modes are not used here.

`dual_agn.provenance.verify_imported_scorer` verifies the actual imported file,
its Git checkout, the full pinned revision, and a clean working tree. A
configured folder at the right revision cannot conceal an import from another
environment. Scoring entry points, calibration and control workers perform the
check. Dry-run plans and discovery-only runs do not score sources.

The existing calibration JSON records the verified scorer revision for new
runs. Existing historical outputs are not retroactively assigned provenance
that they did not record.
