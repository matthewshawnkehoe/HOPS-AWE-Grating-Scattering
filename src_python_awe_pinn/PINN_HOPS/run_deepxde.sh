#!/bin/sh
# DeepXDE runs: dielectric (tanh, DeepXDE's adaptive "LAAF-10 tanh"), gold (tanh)
for spec in "dielectric|tanh" "dielectric|LAAF-10 tanh" "gold|tanh"; do
  c=${spec%%|*}; a=${spec#*|}
  f="results/variants/deepxde_${c}_$(echo $a | tr ' ' '_').json"
  [ -f "$f" ] || python -u deepxde_solver.py --case $c --activation "$a"
done
