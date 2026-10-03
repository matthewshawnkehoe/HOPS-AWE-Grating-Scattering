#!/bin/sh
# all least-squares PINN results (finished parts are skipped); ~1 h on 2 busy cores
export OPENBLAS_NUM_THREADS=1
[ -f results/lsq/conv.json ] || python -u compare_lsq.py --parts conv
[ -f results/lsq/coords.json ] || python -u compare_lsq.py --parts coords
[ -f results/lsq/map_dielectric_q1.json ] || python -u compare_lsq.py --parts map --map-scenarios dielectric
[ -f results/lsq/map_gold_q1.json ] || python -u compare_lsq.py --parts map --map-scenarios gold
[ -f results/inverse/inverse.json ] || python -u inverse_demo.py
[ -f results/lsq/points.json ] || python -u compare_lsq.py --parts points
