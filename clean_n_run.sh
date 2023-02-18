#!/bin/bash

# clean folder
for filename in "*.png" "*.vtu" "*.pvd" "*.pvtu" "cost_history.txt" "./shots/*" "shot_number_level_set_*" "./results/*" 
do
    rm -r $filename
done
# run experiment

# for code in "./run_forward_sharp_interface_generic.py" "./run_sharp_interface_optimization_generic.py"
# do
#     mpiexec -n $n_cores python $code
# done


