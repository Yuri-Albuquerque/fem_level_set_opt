#!/bin/sh

export PYTHONPATH=${PYTHONPATH}:/home/yuri/project:/home/yuri/project/fem_level_set_opt
export FIREDRAKE_CACHE_DIR=/home/yuri/project/fem_level_set_opt/cache_dir/firedrake
export PYOP2_CACHE_DIR=/home/yuri/project/fem_level_set_opt/cache_dir/pyop2
export FIREDRAKE_TSFC_KERNEL_CACHE_DIR=/home/yuri/project/fem_level_set_opt/cache_dir/tfsc
export OMP_NUM_THREADS=1 
