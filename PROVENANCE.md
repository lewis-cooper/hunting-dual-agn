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

Imaging validates exactly two prepared components, distinct Tractor identities,
unique primary/secondary coordinate matches and consistent adopted system
redshift metadata before scoring. Both components must explicitly satisfy
reference-band S/N >= 10, |b| >= 25 degrees and PSF morphology in the selection
mask itself. The cutflow uses those same quality masks.

Discovery still analyses the two brightest r-band PSF sources in an aperture.
`n_psf_components_within_aperture` records the full PSF source count, including
additional sources, in the imaging audit and combined final catalogue. It is
missing for Gaia-only systems whose apertures were not evaluated by imaging;
a missing count must not be read as two. Additional sources are not an automatic
exclusion. This remains a two-component photometric analysis of such fields.

Scoring writes its verified import record to `scorer_provenance.json` in the
working directory. Spawned blend workers verify their own imports too.

Confirmed lenses are retained and flagged. Removal requires a verified pair
identity and a confirmed non-lens contaminant classification. Uncertain and
conflicting cases remain flagged. The reference CSVs and matching rules are
documented in [data/literature/README.md](data/literature/README.md).

Visual inspection remains outside the automatic workflow. Candidate counts and
the calibrated threshold are measured from supplied data and are not fixed to
the historical reference values.

## Published output snapshot

The tracked `products/` and `figures/` contain the completed 26 November 2026
candidate workflow from supplied survey catalogues through discovery, component
preparation, scoring, calibration, morphology, literature auditing and plots.
Discovery recovers 712 Gaia-neighbour entries; imaging passes 1,600 systems to
component scoring. The calibration gives ln R >= -6.81 from 109 systems using
the verified pinned scorer.

The method selects 442 unique systems before literature exclusions. Removing
five confirmed non-lens contaminants leaves 437 unique candidates: 98 resolved,
39 blended and 335 imaging selections, with 35 shared between the searches.
After the 27 November literature-coordinate refresh, nine candidates match
confirmed dual AGN and 26 are flagged as confirmed lenses. The refresh adds
independent component identifiers for J0402-3237 and target coordinates for
J0400+2754 and J2332+2243. The confirmed-dual comparison counts a published
position within 1 arcsec of either candidate component; this also accepts
J1552+3009, whose second absolute coordinate is unavailable. Contaminant
removal still requires the existing pair evidence. Candidate membership and
scores are unchanged. Visual-inspection exclusions are not applied. The six products and
18 PDF/PNG figure files were regenerated; the final workflow verification
passed all checks, including component identities, quality cuts, multiplicity
metadata, redshift consistency and literature exclusions.

The independent chance-alignment control experiment is not repeated by this
candidate run. Its rates must be distinguished from the final sample counts.

Regenerated products and figures replace these local files and can be committed
together once validated. Analysis scripts do not automatically commit or upload
data. Native inputs, cutouts, working caches and credentials remain excluded.
