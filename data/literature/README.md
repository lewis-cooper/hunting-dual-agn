# Literature reference snapshots

These CSVs are curated publication data used by the offline crossmatch stage.
The snapshot was assembled on 23 November 2026 from the saved confirmed-dual
ledger and published target/classification tables. It contains no project
candidate identifiers, private survey inputs or candidate-derived coordinates.

`known_binaries.csv` contains 55 distinct confirmed dual-AGN systems within
0 <= redshift <= 4 and 0 < angular separation <= 1.5 arcsec. A compilation's
explicit confirmed label is accepted. Disputed dual/lens cases are excluded
from this confirmed set. J0400+2754 and J2332+2243 have no usable precise
coordinates and remain present with blank positions. Coordinates are never
invented from an abbreviated J name.

`known_lenses_contaminants.csv` contains 1,617 published pair/classification
records. Its sources include the GraL X Table A.1 lens images, Lemon et al.
(2023), SQLS follow-up outcomes, GMP/ERIS spectroscopy, Shen et al. (2023)
and the earlier position/classification audit. Every record cites its source.
This is a curated comparison catalogue, not an exhaustive census. Several
records can refer to the same system: lenses with several images have distinct
image-pair rows, and papers can repeat or qualify classifications. Counts must
deduplicate matched candidate IDs rather than count reference rows. Contaminant
references are not restricted by the dual-reference redshift/separation window;
a close pair can comprise two images of a wider multiple lens.

## Columns

| Columns | Meaning |
| --- | --- |
| `reference_id`, `system_id`, `system_name` | Unique reference-record key, literature system/alias label and published name. Different aliases can describe the same physical system. |
| `ra_deg`, `dec_deg` | Published component or system position, in J2000 degrees. Its role and origin are given by `astrometry_source`. |
| `secondary_ra_deg`, `secondary_dec_deg` | Second component position when published or derivable from a published component position, separation and position angle. Derived positions are explicitly labelled. |
| `redshift`, `secondary_redshift` | Published system/component redshifts where available. Blank means unavailable; these values do not gate matches. |
| `separation_arcsec`, `position_angle_deg` | Published or explicitly derived pair geometry; PA runs north through east. |
| `primary_gaia_source_id`, `secondary_gaia_source_id` | Exact Gaia DR3 identifiers when supplied. Read them as strings. |
| `classification` | Adopted classification, retaining uncertain and disputed states. |
| `author_year`, `publication_url`, `source_table` | Publication attribution and the location of the measurements/classification. Multiple confirming sources are retained for duals. |
| `astrometry_source`, `separation_source` | Origin of the positions and separation, including derived values. |
| `match_rule`, `position_tolerance_arcsec`, `separation_tolerance_arcsec` | Documented pair-identification rule and tolerances. |
| `notes` | Relevant qualifications, conflicts and evidence. |

## Pair identification

The default 1 arcsec radius is only a finder. The matcher accepts identical
unordered pairs of two distinct Gaia IDs or both component positions within
the listed tolerance, allowing exchanged primary/secondary labels. Gaia and
precise published component pairs normally use 0.1 arcsec. ERIS positions
derived from measured separation/PA use 0.2 arcsec, accommodating the
documented J0008-3655 centroid/separation difference.

`published_pair_geometry` requires a component/system position and pair
separation to agree within the listed tolerances; when PA is available, the
unordered axes must agree within 5 degrees. It applies to documented
spectroscopic pair identifications with no separate absolute component
positions. SQLS quasar–star positions use up to 0.4 arcsec because the SDSS
positions and Gaia components differ, and also require separation agreement
within 0.15 arcsec. The reviewed J0216+1430 and J2335+3201 dual identities use
0.1 arcsec in position and 0.02 arcsec in separation.

`reviewed_system_geometry` also permits a published system position to agree
with the pair midpoint. This is used for the previously reviewed J2302-4154
projection: its rounded system position agrees within 0.3 arcsec and its
reported 1.22 arcsec separation must agree within 0.15 arcsec. Its published
foreground redshift depends on a single-line identification; that qualification
remains in the notes. `position_only` records remain possible associations
unless both published Gaia IDs establish the pair.

## Classification policy

Only verified `confirmed_lens`, `confirmed_quasar_star`,
`confirmed_quasar_galaxy`, `confirmed_stellar_pair`, and
`confirmed_projected_quasar_pair` matches are removed. Likely lenses,
lens candidates, possible quasar–star pairs and unclassified quasar pairs are
retained and flagged. Conflicting confirmed dual/contaminant classifications
require review and do not trigger automatic removal.

J0823+2418 is retained as `likely_lens`, following the qualified dedicated
study rather than promoting a catalogue lens label to an unambiguous
confirmation. J0749+2255 and J0841+4825 retain disputed dual/lens states.
J1059+0622 retains the qualified lens-candidate interpretation in the source
discussion. Shen's `QS_PM` and `QS_PCA` labels are indirect stellar flags,
not automatic confirmations of stellar contamination.

The GraL image catalogue is attributed to
[Ducourant et al. (2026)](https://arxiv.org/abs/2603.12810),
[VizieR J/A+A/707/A345, Table A.1](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/707/A345).
Only individually listed quasar-image components are paired; deflectors are
excluded. Entries without separate image positions retain only the published
system position and cannot by themselves trigger an automatic exclusion.
Publication URLs and table locations in each CSV provide the remaining
record-level provenance. Updating these curated inputs can change the output
counts; the run summary records their SHA-256 hashes.
