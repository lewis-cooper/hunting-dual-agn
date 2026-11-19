#!/bin/sh
set -eu
paper_python="${PAPER1_PYTHON:-python}"
"$paper_python" -m unittest discover -s verification -p 'test_*.py'
"$paper_python" verification/verify_redshift_uncertainty.py
"$paper_python" verification/verify_input_contracts.py
for stage in \
    prepare_one_arcsecond_inputs score_one_arcsecond_resolved \
    calibrate_resolved_threshold score_one_arcsecond_unresolved \
    select_one_arcsecond_candidates select_imaging_precolour \
    prepare_imaging_components score_imaging_components \
    select_imaging_candidates \
    combine_candidate_catalogues prepare_photometry_inventory \
    plot_one_arcsecond_results plot_one_arcsecond_unresolved \
    plot_imaging_results plot_joint_results
 do
    "$paper_python" -u -m "scripts.$stage"
 done
"$paper_python" verification/verify_example_run.py
