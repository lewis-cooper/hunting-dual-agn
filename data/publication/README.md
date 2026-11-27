# Publication decisions

`exclusions.csv` records the corrupted-data artefact removed after the automatic
selection. Both Gaia component identifiers must agree before an exclusion is
applied. An ordinal candidate ID alone is not sufficient. This curation is
separate from the automatic scoring and cutflow.

`euclid_visual_classifications.csv` preserves the existing author-confirmed
image review: the usable systems show two sources except J094112+621719, which
is a potential blend. These are visual multiplicity assessments, not
confirmations of two AGN. Footprint membership and usable-image counts are
recalculated independently by `recount_euclid.py` from supplied local MOCs and
the cached VIS pixels. A missing image cache or review is an error, not an
assumed non-detection.

The image cache and footprint FITS files are not included here. Their paths are
supplied through the optional Euclid configuration when exporting the catalogue.
