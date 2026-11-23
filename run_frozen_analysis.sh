#!/bin/sh
set -eu
paper_python="${PAPER1_PYTHON:-python}"
exec "$paper_python" -u -m scripts.run_analysis --stage vetting "$@"
