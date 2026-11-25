# Validation

Run the tests from the repository root after installing the pinned scorer:

```sh
python -m unittest discover -s tests -p 'test_*.py'
python tests/verify_redshift_uncertainty.py
python tests/verify_input_contracts.py --models-only
```

The tests exercise discovery, component association, native photometry, scoring
selection rules, scorer provenance, literature matching and workflow ordering.
They use small fixtures and mocked services; passing tests does not establish
the number of objects in a real survey run.

The original survey exports have been checked against the historical discovery
identities: 1,269 qualifying associations become 865 catalogue entries and 712
position groups. A live neighbour-count lookup also recovered those groups.
These are initial inputs to vetting, not a required count or final dual-AGN
classification. Different supplied populations can produce different counts.

After a complete candidate run, check its catalogue and audit consistency:

```sh
python tests/verify_example_run.py --pipeline both
```

The final CSV, score audits, cut tables, literature associations and calibration
must agree. Scorer provenance is recorded in the calibration JSON for new runs.
Historical reference files did not record all imported-package and working-tree
details and are not retroactively certified by the current environment check.

Visual inspection, Euclid image classifications and interpretation of the
control rates require separate review. Control acceptance probabilities are not
measurements of final-sample contamination or purity. The geometric field-density
diagnostic does not supply a measured close-pair detection completeness.
