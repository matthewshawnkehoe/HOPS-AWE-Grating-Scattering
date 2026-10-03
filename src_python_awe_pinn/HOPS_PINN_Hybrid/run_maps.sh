#!/bin/sh
# all map comparisons (finished scenarios are skipped; truth maps are cached)
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -u compare_maps.py "$@"
