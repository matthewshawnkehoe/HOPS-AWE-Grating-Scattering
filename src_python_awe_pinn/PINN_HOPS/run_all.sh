#!/bin/bash
# Full comparison campaign (several hours on 2 CPU cores). The run_*.sh scripts skip finished runs.
cd "$(dirname "$0")"
# 1. first study: gradient-trained PINN
python compare_point.py --budget standard
python compare_refl_map.py --scenario dielectric --budget standard
python compare_refl_map.py --scenario gold --budget standard
# 2. second study: better PINN variants
./run_variants.sh            # activations (tanh/sin/LAAF), I-PINN, LSQ output layer, LSGD  (quick budget)
./run_deepxde.sh             # the same problem in DeepXDE
./run_lsq.sh                 # least-squares interface PINN: convergence, coordinates, maps, inverse problem, points
python timing_benchmark.py   # cost per step, single thread
python summarize_variants.py # results/variants/summary.md, summary.png
