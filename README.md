# Hunting Dual AGN

Analysis code for *Hunting arcsecond-scale dual and lensed quasars in wide-field
surveys*, by Lewis M. Cooper, Vasily Belokurov and Sergey E. Koposov.

The repository searches for close quasar pairs through Gaia neighbours and
Legacy Surveys imaging. It builds component measurements, scores resolved
companions and blended photometry, calibrates the selection threshold, measures
morphology, and combines the candidates with a literature audit. Confirmed
lenses are retained and flagged; confirmed non-lens contaminants are removed.

Users supply the quasar catalogue, survey extracts, coverage records and FITS
images. Observations, model assets and credentials are not distributed here.
Generated products and figures are tracked. The automatic full pipeline run
contains 437 candidates, retaining confirmed lenses and leaving visual inspection
outside the automatic selection; see [provenance](PROVENANCE.md). The
[publication catalogue](catalogue/README.md) contains 436 systems after the
explicit corrupted-data exclusion: 201 main and 235 additional. It is supplied
as a combined CSV and two disjoint tier tables, with uniform coordinate names
and preserved original identifiers. The `tier`
column labels systems with two distinct Gaia detections as `main` and one-Gaia
imaging systems as `additional`. Tier assignment does not change selection.
The optional chance-alignment experiment has its own control-data acquisition.

Start with `run_analysis.py --pipeline both`; either `gaia_neighbour` or `imaging`
can also be selected. Input configuration and installation are required first.

- [Data requirements and catalogue definitions](DATA.md)
- [Installation, configuration and run commands](SETUP.md)
- [Pinned scorer and scientific selection rules](PROVENANCE.md)
- [Tests and validation limits](VALIDATION.md)
- [Chance-alignment controls](chance_alignments/README.md)
- [Publication catalogue export](SETUP.md#export-the-publication-catalogue)

Please cite the accompanying paper when using this analysis.
