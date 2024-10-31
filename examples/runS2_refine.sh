#!/bin/bash

exec=/Users/jinyunlin/QFE-Research/QFEgeometry/bin/S2_refine
q=4
N_iter=1000
stepsize=1
data_dir=/Users/jinyunlin/QFE-Research/QFEgeometry/data/experiment/q${q}
#for L in 8 10 24 36 52 72 84
#for L in 4 8 10 16 20 24 32 36 48 52 64 72 84 96
for L in 84
do 
orbit_dir=/Users/jinyunlin/QFE-Research/QFEgeometry/data/q${q}/xivec_${L}.dat
mkdir -p ${data_dir}
${exec} ${L} ${N_iter} ${stepsize} ${data_dir} ${q} ${orbit_dir} >> ${data_dir}/L${L}.out
done
