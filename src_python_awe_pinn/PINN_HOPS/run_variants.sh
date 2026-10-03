#!/bin/sh
# gradient-trained PINN variants (activation study, I-PINN, LSQ output layer, LSGD); finished runs are skipped
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for v in tanh sin laaf ipinn_tanh_sin tanh+lsq lsgd; do
  for c in dielectric gold; do
    [ -f "results/variants/act_${c}_${v}.json" ] || python study_activations.py --cases $c --variants $v
  done
done
