# Scorer and selection provenance

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

## Selection and literature policy

Either model may select a resolved companion; imaging requires both components
to pass within the same model variant. Threshold calibration precedes the final
secondary-PM cut and keeps its calibration-only upper limit of 70. Two-Gaia
systems require secondary PM significance below 3 when finite; missing PM is
retained and flagged. Imaging also requires finite primary PM significance below
3. The DCR angle threshold of 20 degrees applies once to one-Gaia imaging only.
Morphology uses R_perp >= 1.05. Quaia redshift errors are propagated in scoring.

Confirmed lenses are retained and flagged. Removal requires a verified pair
identity and a confirmed non-lens contaminant classification. Uncertain and
conflicting cases remain flagged. The reference CSVs and matching rules are
documented in [data/literature/README.md](data/literature/README.md).

Visual inspection remains outside the automatic workflow. Candidate counts and
the calibrated threshold are measured from supplied data and are not fixed to
the historical reference values.

## Published output snapshot

The tracked `products/` and `figures/` contain the completed 24 November 2026
reference results: 411 unique candidates, score audits for 712 Gaia-neighbour
inputs and 1,600 imaging systems, two cutflows, the threshold calibration, and
18 PDF/PNG figure files. The calibration gives ln R >= -6.81 from 109 systems.

This snapshot applied saved visual-inspection exclusions and removed confirmed
lenses as well as other contaminants. It therefore predates the current policy
above, which retains confirmed lenses and leaves visual inspection outside the
automatic workflow. These files are reference results, not a completed rerun
of the current workflow. The original calibration JSON is preserved without
retrospectively adding an imported-scorer verification record.

Regenerated products and figures replace these local files and can be committed
together once validated. Analysis scripts do not automatically commit or upload
data. Native inputs, cutouts, working caches and credentials remain excluded.
